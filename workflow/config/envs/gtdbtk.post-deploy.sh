root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)
config_path="${root_path}"/workflow/config
db_root=$(awk -F ': ' '/db_root:/{print $2}' "${config_path}"/config.yaml)
gtdbtk_db=$(awk -F ': ' '/    gtdbtk:/{print $2}' "${config_path}"/config.yaml)

  {
    echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (gtdbtk)"
    echo "---------------------------"
  } >> "${root_path}"/logs/env.log

# Check whether the gtdbtk database already exists
if [ -d "${db_root}"/"${gtdbtk_db}" ];then
  conda env config vars set GTDBTK_DATA_PATH="${db_root}/${gtdbtk_db}"
  echo "GTDB-Tk database already exists, location：${db_root}/${gtdbtk_db}" >> "${root_path}"/logs/env.log
  echo "To check the database completeness, run the following command: gtdbtk check_install" >> "${root_path}"/logs/env.log
else
  echo "The GTDB-Tk database doesn't exist. Use the following command to build it yourself:" >> "${root_path}"/logs/env.log
 ### 注意，EOF语句不允许前面有tab或空格，所以此处不缩进
cat >> "${root_path}"/logs/env.log <<EOF
    GTDB-Tk v2.3.2 requires ~78G of external data which needs to be downloaded and extracted.
    This can be done automatically, or manually.

    Automatic:
        1. Run the command "download-db.sh" to automatically download and extract to:
            {$CONDA_PREFIX}/share/gtdbtk-2.3.2/db/
    Manual:
        1. Manually download the latest reference data:
            wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
        2. Extract the archive to a target directory:
            tar -xvzf gtdbtk_r214_data.tar.gz -C "${db_root}/${gtdbtk_db}" --strip 1 > /dev/null
            rm gtdbtk_r214_data.tar.gz
        3. Set the GTDBTK_DATA_PATH environment variable by running:
            conda env config vars set GTDBTK_DATA_PATH=${db_root}/${gtdbtk_db}
EOF
exit 1
fi
