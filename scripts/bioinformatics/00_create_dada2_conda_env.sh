module load python/miniforge-25.3.0

conda create -n amplicon -c conda-forge -c bioconda cutadapt

conda activate amplicon
conda init bash
source ~/.bashrc

conda activate amplicon
