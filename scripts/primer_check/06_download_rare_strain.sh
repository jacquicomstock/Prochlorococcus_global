esearch -db biosample -query "SAMN38315731" | \
elink -target nuccore | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element AccessionVersion Title

esearch -db biosample -query "SAMN38315732" | \
elink -target nuccore | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element AccessionVersion Title

esearch -db biosample -query "SAMN38315733" | \
elink -target nuccore | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element AccessionVersion Title

esearch -db biosample -query "SAMN38315735" | \
elink -target nuccore | \
efetch -format docsum | \
xtract -pattern DocumentSummary -element AccessionVersion Title

efetch -db nucleotide -id CP139303.1 -format gbwithparts > rare_strain_genbank/MIT1223.gb
efetch -db nucleotide -id CP139302.1 -format gbwithparts > rare_strain_genbank/MIT1300.gb
efetch -db nucleotide -id CP139301.1 -format gbwithparts > rare_strain_genbank/MIT1307.gb
efetch -db nucleotide -id CP139304.1 -format gbwithparts > rare_strain_genbank/MIT1341.gb
