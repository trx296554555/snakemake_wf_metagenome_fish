root_path=$(dirname `dirname \`dirname "$CONDA_PREFIX"_\``)
env_name=$(basename $CONDA_PREFIX)

# https://www.nature.com/articles/s41592-022-01419-0
# Critical Assessment of Metagenome Interpretation: the second round of challenges
# https://www.nature.com/articles/s41592-022-01431-4
# In CAMI2 Genome binning challenge metabinner showed good performance
# A stand-alone ensemble binning method

ln -s $CONDA_PREFIX/bin/scripts/* $CONDA_PREFIX/bin/

{
  echo "Make additional adjustments for the post-deployment of the Conda environment ${env_name} (metabinner)"
  echo "Binning the assembled Contigs using metabinner."
  echo "---------------------------"
} >> "${root_path}"/logs/env.log