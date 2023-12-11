import os, subprocess
import asyncio
from async_download import main_download


def download_ref_fa_file(ref_dict, db_dir):
    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    else:
        already_downloaded = os.listdir(db_dir)
        need_download = {}
        for species in ref_dict:
            fna = f'{str(species).replace(" ", "_")}.fna.gz'
            md5 = f'{str(species).replace(" ", "_")}.fna.gz.md5'
            if fna in already_downloaded:
                pass
            else:
                need_download[os.path.join(db_dir, fna)] = ref_dict[species]
            if md5 in already_downloaded:
                pass
            else:
                need_download[os.path.join(db_dir, md5)] = os.path.dirname(ref_dict[species])+'/md5checksums.txt'
    # print(need_download)
    asyncio.run(main_download(need_download))



def ref_fa_file_md5_check(ref_dict, db_dir):
    try:
        with open(os.path.join(db_dir, 'md5_check.log'), 'r') as f:
            for line in f:
                line = line.split(':')
                if 'OK' in line[1] or '成功' in line[1]:
                    ref_dict.pop(line[0].split('.')[0], None)
    except FileNotFoundError:
        pass
    if not ref_dict:
        print('All ref files are downloaded successfully!')
        return None
    for species in ref_dict:
        fna = f'{str(species).replace(" ", "_")}.fna.gz'
        md5 = f'{str(species).replace(" ", "_")}.fna.gz.md5'
        raw_file = os.path.basename(ref_dict[species])
        with open(os.path.join(db_dir, md5), 'r') as f:
            for line in f:
                if raw_file in line:
                    raw_md5 = line.split()[0]
                elif fna in line:
                    raw_md5 = line.split()[0]
                else:
                    pass
        with open(os.path.join(db_dir, md5), 'w') as f:
            f.write(f'{raw_md5}  {fna}\n')
        command = f'cd {db_dir} && md5sum -c {os.path.join(db_dir, md5)} >> {os.path.join(db_dir, "md5_check.log")}'
        if subprocess.check_output(command, shell=True, executable="/bin/bash", text=True):
            raise Exception(f'{fna} is not downloaded successfully!')