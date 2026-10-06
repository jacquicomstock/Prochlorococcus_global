# The purpose of this script is to check if there is variability in the 27F or 338R primer binding sites 
# between Prochlorococcus clades that could potentially bias amplification and sequencing output.

#create new conda env
conda create -n entrez -c bioconda entrez-direct -y
conda activate entrez

# see how many potentially useful Prochlorococcus 16S record NCBI has prior to downloading
esearch -db nucleotide \
  -query '"Prochlorococcus"[Organism] AND ("16S ribosomal RNA"[Title] OR "16S rRNA"[Title]) AND 1200:1700[SLEN]' \
  | xtract -pattern ENTREZ_DIRECT -element Count

  #Download all 2016 sequences
  esearch -db nucleotide \
  -query '"Prochlorococcus"[Organism] AND ("16S ribosomal RNA"[Title] OR "16S rRNA"[Title]) AND 1200:1700[SLEN]' \
  | efetch -format fasta \
  > Prochlorococcus_full16S.fasta

  #extract 27f sequences
  mafft --auto Prochlorococcus_full16S.fasta > Prochlorococcus_full16S_aligned.fasta
