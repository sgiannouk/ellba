### This script is designed to facilitate the download of all essential reference files required for the ELLBA workflow

import os
import subprocess
from datetime import datetime
startTime = datetime.now()



main_dir 	= os.path.dirname(os.path.realpath(__file__))
dbsnp_dir 	= os.path.join(main_dir, "dbSNP")
gatk_dir 	= os.path.join(main_dir, "gatk_subfiles", "known_sites")

### REFERENCE FILES
print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Step1. DOWNLOADING REFERENCE GENOME, TRANSCRIPTOME AND ANNOTATION FILES')

# Gencode Reference Human Genome
ref_genome = os.path.join(main_dir, "GRCh38.primary_assembly.genome.fa")
ref_genome_dict = os.path.join(main_dir, "GRCh38.primary_assembly.genome.dict")
gencode_ref_genome_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz"
subprocess.run(f"wget -O - {gencode_ref_genome_url} | pigz -d -c > {ref_genome}", shell=True, check=True)
subprocess.run(f"samtools faidx {ref_genome}", shell=True, check=True)
# Creating sequence dictionary for the reference genome
subprocess.run(f"picard CreateSequenceDictionary R={ref_genome} O={ref_genome_dict}", shell=True, check=True)


# Gencode v35 Reference Human Genome Annotation
ref_annotation = os.path.join(main_dir, "gencode.v35.primary_assembly.annotation.gtf")
gencode_ref_annotation_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz"
subprocess.run(f"wget -O - {gencode_ref_annotation_url} | pigz -d -c > {ref_annotation}", shell=True, check=True)

# Gencode v35 Reference Human Transcriptome
ref_tran = os.path.join(main_dir, "gencode.v35.transcripts.orig.fa")
ref_transcriptome = os.path.join(main_dir, "gencode.v35.transcripts.fa")
gencode_ref_transcriptome_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz"
subprocess.run(f"wget -O - {gencode_ref_transcriptome_url} | pigz -d -c > {ref_tran}", shell=True, check=True)
# Manipluating the header of the fasta file
maniplulate_command = "sed -n -e \'s/^>\\(ENST[^|]*\\)|.*$/>\\1/p\' -e \'/^>/!p\'"
subprocess.run(f"{maniplulate_command} {ref_tran} > {ref_transcriptome}", shell=True)
os.remove(ref_tran)


### REFERENCE DATABASES
print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Step2. DOWNLOADING REFERENCE DATABASES')

# DBSNP
if not os.path.exists(dbsnp_dir): os.makedirs(dbsnp_dir)

dbsnp = os.path.join(dbsnp_dir, "common_all_20180418.chr.vcf")
dbsnp_comp = os.path.join(dbsnp_dir, "common_all_20180418.chr.vcf.gz")
dbsnp_orig = os.path.join(dbsnp_dir, "common_all_20180418.chr.orig.vcf")
dbsnp_url = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz"
subprocess.run(f"wget -O - {dbsnp_url} | pigz -d > {dbsnp_orig}", shell=True, check=True)

print('You need to wait juust a little bit more to manipulate the dbSNP file...')
# Add chr infront of each chromosome in the CHROM column 
maniplulate_command = "awk \'{if($0 !~ /^#/) print \"chr\"$0; else print $0}\'"
subprocess.run(f"{maniplulate_command} {dbsnp_orig} > {dbsnp}", shell=True)
# Compressing the vcf and indexing it
subprocess.run(f"bgzip {dbsnp} && tabix -p vcf {dbsnp_comp}", shell=True, check=True)

# Converting dbSNP from vcf to tab-delimited file
dbsnp_tab_orig = os.path.join(dbsnp_dir, "common_all_20180418.chr.orig.tab")
dbsnp_tab = os.path.join(dbsnp_dir, "common_all_20180418.chr.tab")
subprocess.run(f"gatk VariantsToTable -V {dbsnp_orig}  -F CHROM -F POS -F ID -F REF -F ALT -O {dbsnp_tab_orig}", shell=True, check=True)
maniplulate_comm = "awk \'NR==1; NR>1 && $0 !~ /^#/{print \"chr\"$0}\'"
subprocess.run(f"{maniplulate_comm} {dbsnp_tab_orig} > {dbsnp_tab}", shell=True)
os.remove(dbsnp_orig)
os.remove(dbsnp_tab_orig)


# GATK
if not os.path.exists(gatk_dir): os.makedirs(gatk_dir)

thousandG_omni = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
thousandG_omni_phase1 = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
known_indels = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
mills_and_thousandG = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
all_urls = [thousandG_omni, thousandG_omni_phase1, known_indels, mills_and_thousandG]

for url in all_urls:
	file_name = os.path.basename(url)
	save_output = os.path.join(gatk_dir, file_name[:-3])
	subprocess.run(f"wget -O - {url} | pigz -d -c > {save_output}", shell=True, check=True)
	subprocess.run(f"gatk IndexFeatureFile -I {save_output}", shell=True, check=True)

print("\n\n\n\n\n\nFinally, all files are now download and you can run the ELLBA workflow :)")