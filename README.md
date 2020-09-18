# dbCOG
Build COG database
## URL of COG database
[https://www.ncbi.nlm.nih.gov/research/cog/](https://www.ncbi.nlm.nih.gov/research/cog/)  
[https://www.ncbi.nlm.nih.gov/research/cog-project/](https://www.ncbi.nlm.nih.gov/research/cog-project/)  
[ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/)  

## Database directory
[Readme.2020-09-15.txt](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/Readme.2020-09-15.txt)  
[cog-20.cog.csv](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv)  
[cog-20.def.tab](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab)  
[cog-20.org.csv](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.org.csv)  
[cog-20.tax.csv](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.tax.csv)  
[fun-20.tab](ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab)  

## Building Requirements
[biopython](https://biopython.org/)

## Download database
<pre><code>
wget -c ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
git clone https://github.com/zxgsy520/dbCOG.git
python download_cog.py cog-20.cog.csv >COG.pep.fa
#of
python download_cog_seqs.py cog-20.cog.csv --thread>COG.pep.fa #Download data using multiple processes
<pre><code>
