#create new conda env
conda create -n phyloassigner \
  -c conda-forge \
  -c bioconda \
  python=3.9 \
  pandas \
  biopython \
  hmmer \
  pplacer \
  gappa \
  taxtastic \
  raxml \
  -y

  #activate
  conda activate phyloassigner
