import os
import subprocess

sample_dir = os.path.dirname(os.path.realpath(__file__))
clin_data = os.path.join(sample_dir, "clinical_data.tsv")


# Stavros Giannoukakos
""" Omitting fastq.sra_bd (bioinfokit.analys toolkit) 
and downloading from amazon AWS server"""

#Version of the program
__version__ = "0.0.3"


import glob, os
import subprocess
from datetime import datetime
startTime = datetime.now()



sra_clindata  = f"{sample_dir}/SraRunTable.txt"
sra_list 	  = f"{sample_dir}/sra_accession.txt"
stats_dir 	  = f"{sample_dir}/stats"
threads 	  = 2


def get_sras():
	print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} Downloading the requested fastq files: in progress..')

	downloaded_files = [os.path.basename(files).split(".")[0] for files in glob.glob(os.path.join(sample_dir, "*.fastq"))]
	total = sum(1 for line in open(sra_list))

	with open(sra_list) as fin:
		for i, line in enumerate(fin, 1):
			srr = line.strip()
			if srr not in downloaded_files:	
				print(f"{i}/{total}  Downloading: {srr}")
				subprocess.run(f'wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/{srr}/{srr} -O {sample_dir}/{srr}.sra', shell=True)
				subprocess.run(f'fastq-dump --split-files {sample_dir}/{srr}.sra', shell=True)
				subprocess.run(f'rm {sample_dir}/{srr}.sra', shell=True)
	return

def obtain_sra_stats():
	print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} Obtaining stats of the downloaded fastq files: in progress..')
	#### Check number of bases
	if not os.path.exists(stats_dir): os.makedirs(stats_dir)

	### Existing stat files
	stat_files = [os.path.basename(files).replace(".txt","") for files in glob.glob(f"{stats_dir}/*.txt")] 

	for i, files in enumerate(glob.glob(f"{sample_dir}/*.fastq")):
		file_name = files.split("/")[-1].split(".")[0]
		if file_name not in stat_files:
			os.system(f'reformat.sh in={files} 2>> {stats_dir}/{file_name}.txt')
	return

def verify_sras_n_find_missing_data():
	print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} Verifying missing or corrupted fastq files: in progress..')
	### Obtain all SRA files from sra accession list
	total_sras = []
	with open(sra_list) as fin:
		for line in fin:
			total_sras.append(line.strip())
	print("Total SRA files:",len(total_sras))


	### Obtain bp info from clinical data file
	clinical_data = {}
	with open(sra_clindata) as cin:
		for line in cin:
			if not line.startswith("Run"):
				sample = line.strip().split(",")[0]
				bases = int(line.strip().split(",")[3])
				clinical_data[sample] = bases
	

	downloaded_files = {}
	for stat_files in glob.glob(f"{stats_dir}/SRR*_1.txt"):
		with open(stat_files) as sin:
			sample_name = os.path.basename(stat_files).split(".")[0].split("_")[0]
			for line in sin:
				if line.startswith("Input:"):
					downloaded_files[sample_name] = int(line.strip().split()[-2])
		with open(stat_files.replace("_1.", "_2.")) as sin2:
			sample_name = os.path.basename(stat_files).split(".")[0].split("_")[0]
			for line2 in sin2:
				if line2.startswith("Input:"):
					downloaded_files[sample_name] += int(line2.strip().split()[-2])

	i = 0
	for sample_name, bases in downloaded_files.items():
		if sample_name in clinical_data:
			# print(i, sample_name, f"reference:{clinical_data[sample_name]}", f"fastq has {bases}")
			if bases != clinical_data[sample_name]:
				i+=1
				print(i, sample_name, f"reference:{clinical_data[sample_name]}", f"fastq has {bases}")

	print()
	print()
	print("\tIn total", len(list(set(total_sras) - set(downloaded_files))), "samples weren't downloaded")
	for j, samps in enumerate(list(set(total_sras) - set(downloaded_files)), 1):
		print("\n\n", j, "Try to download sample:", samps)
		subprocess.run(f'prefetch {samps}',shell=True)
		subprocess.run(f'vdb-validate {samps}',shell=True)
		subprocess.run(f'fasterq-dump --threads {threads} --progress {samps}',shell=True)
	return

def gz_fastq():
	print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} Compressing all fastq files: in progress..')

	fastq_files = ' '.join(glob.glob(f'{sample_dir}/*.fastq'))
	subprocess.run(f'pigz --processes {threads} {fastq_files}', shell=True)
	return

def rename():
	print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")} Renaming all fastq.gz files: in progress..')

	to_rename = {}
	with open(clin_data) as sin:
		for line in sin:
			if not line.startswith("SampleName"):
				run = line.strip().split()[2]
				sampleName = line.strip().split()[0]
				to_rename[run] = sampleName

	for file in os.listdir(sample_dir):
		if file.endswith(".fastq.gz"):
			sample_path = sample_dir
			sample_id = file.split("_")[0]
			strand = ["R1" if file.split("_")[1].split(".")[0] == "1" else "R2"][0]
			new_sample_id = to_rename[sample_id]
			old_sample = f"{sample_path}/{file}"
			new_sample = f"{sample_path}/{new_sample_id}.{strand}.fastq.gz"
			subprocess.run(f"mv {old_sample} {new_sample}", shell=True)
	return

def main():

	print("Running SRA toolkit to download the HCC dataset")
	get_sras()
	obtain_sra_stats()
	verify_sras_n_find_missing_data()
	gz_fastq()
	rename()

	print(f'The pipeline finished after {datetime.now() - startTime}')

if __name__ == "__main__": main()
