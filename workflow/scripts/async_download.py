import asyncio
import httpx
from tqdm.asyncio import tqdm_asyncio, tqdm


async def download_file(sem, url, file_path):
    async with sem:
        async with httpx.AsyncClient(timeout=60) as client:
            async with client.stream('GET', url) as response:
                if response.status_code == 200:
                    total_size = int(response.headers.get("Content-Length", 0))
                    with open(file_path, 'wb') as file:
                        async for chunk in tqdm(iterable=response.aiter_bytes(1), desc=file_path, unit='B',
                                                unit_scale=True, unit_divisor=1024, total=total_size):
                            file.write(chunk)
                else:
                    print(f"Failed to download: {url}")


async def main_download(url_dict, max_tasks=10):
    tasks = []
    semaphore = asyncio.Semaphore(max_tasks)
    for file in url_dict:
        tasks.append(download_file(semaphore, url_dict[file], file))
    # await tqdm_asyncio.gather(*tasks, desc='Download Files')
    await asyncio.gather(*tasks)
