#!/bin/bash
#SBATCH --job-name=Download_genomes     
#SBATCH --output=99_logs/03_download_genomes_%j.log                        
#SBATCH --cpus-per-task=2            
#SBATCH --mem=8G                    
#SBATCH --partition=normal
#SBATCH --nodelist=node06

# Script to downlaoad reference genomes as well as potential contaminats (endosymbionts etc)

# Paths
DATA_DIR="/data/onilee/capstone/01_data"

# move into data directory
mkdir -p "$DATA_DIR/03_references"
cd "$DATA_DIR/03_references"

# 1. Download MEAN1 refernece genomes and gff
wget http://www.whiteflygenomics.org/ftp/MEAM1/v1.2/MEAM1_scaffold_v1.2.fa.gz # scafold
wget http://www.whiteflygenomics.org/ftp/MEAM1/v1.2/MEAM1_v1.2.gff3.gz # annotation

# 2. Download mitochondrial genome of MEAN1
wget -O - \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006279.1&rettype=fasta&retmode=text" > NC_006279.1.fa


# 2. Download MEAN1 endosymbiont genomes
wget http://www.whiteflygenomics.org/ftp/MEAM1/endosymbionts/Portiera_meam1_genome_v2.0.fa # portiera
wget http://www.whiteflygenomics.org/ftp/MEAM1/endosymbionts/Hamiltonella_meam1_genome_v2_0.fa # hamiltonella
wget http://www.whiteflygenomics.org/ftp/MEAM1/endosymbionts/Rickettsia_meam1_genome_v2.0.fa # rickettsia

# 3. Merge all contaminnants into one file for easier downstream filtering
cat NC_006279.1.fa Portiera_meam1_genome_v2.0.fa Hamiltonella_meam1_genome_v2_0.fa Rickettsia_meam1_genome_v2.0.fa > MEAM1_contiminants.fa
