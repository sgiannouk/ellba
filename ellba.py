### ANALYSING LIQUID BIOPSY BASED RNASEQ DATA ###
## ESR11 - Stavros Giannoukakos 
## Script: ellba.py

#Version of the program
__version__ = "1.0.1"

# Importing general libraries
import os
import csv
import glob
import math
import json
import random
import shutil
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse, subprocess
import multiprocessing as mp
from datetime import datetime
from collections import Counter
from operator import itemgetter
from multiprocessing import Pool
startTime = datetime.now()



### PATHS
ellba_dir			   = os.path.dirname(os.path.realpath(__file__))
parent_dir			   = os.path.dirname(ellba_dir)
software_dir		   = os.path.join(parent_dir, "software")
references_dir		   = os.path.join(parent_dir, "references")
inputdata_dir		   = os.path.join(parent_dir, "input_data")
rscripts 			   = os.path.join(ellba_dir, "additional_files")
### ADDITIONAL FILES
phix_sequences 		   = os.path.join(rscripts, "phix.fasta")
adapter_sequences 	   = os.path.join(rscripts, "adapter_sequences.fasta")
### FUSION DETECTION
arriba_software		   = os.path.join(software_dir, "arriba_v2.1.0")
arriba 				   = os.path.join(arriba_software, ".", "arriba")
arriba_visualisation   = os.path.join(arriba_software, ".", "draw_fusions.R")
cytobandsGRCh38 	   = os.path.join(arriba_software, "database", "cytobands_hg38_GRCh38_v2.1.0.tsv")
blacklistGRCh38 	   = os.path.join(arriba_software, "database", "blacklist_hg38_GRCh38_v2.1.0.tsv.gz")
knownfusionsGRCh38 	   = os.path.join(arriba_software, "database", "known_fusions_hg38_GRCh38_v2.1.0.tsv.gz")
proteindomainsGRCh38   = os.path.join(arriba_software, "database", "protein_domains_hg38_GRCh38_v2.1.0.gff3")
### ISOFORM DETECTION
ref_salmon 			   = os.path.join(references_dir, "gencode.v35.gffread.fa")
### GATK FILES
gatk_known_sites 	   = os.path.join(references_dir, "gatk_subfiles", "known_sites")
### RNA EDITING
reditools_path		   = os.path.join(software_dir, "REDItools2", "src", "cineca", "reditools.py")
reditools 			   = f"python3 {reditools_path}"
common_dbsnp 	 	   = os.path.join(references_dir, "dbSNP", "common_all_20180418.chr.vcf.gz")
common_dbsnp_tab 	   = os.path.join(references_dir, "dbSNP", "common_all_20180418.chr.tab")
### REFERENCE GENOME and ANNOTATION
referencegenome 	   = os.path.join(references_dir, "GRCh38.primary_assembly.genome.fa")
referenceannot 	   	   = os.path.join(references_dir, "gencode.v35.primary_assembly.annotation.gtf")
### INPUT DATA
input_raw_data		   = os.path.join(inputdata_dir, "raw_data")
input_clinical_data	   = os.path.join(inputdata_dir, "raw_data", "clinical_data.txt")






usage       = "ellba.py [options]"
epilog      = " -- October 2020 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Directory with the raw fastq files
parser.add_argument('-i', '--inputdata', dest='raw_data', default=input_raw_data,
                	help="Path of the directory where the raw fastq or\nfastq.gz data are residing (e.g /home/user/raw_data)")
# Path of the .txt file containing the clinical data
parser.add_argument('-c', '--clindata', dest='clinical_data', default=input_clinical_data,
                	help="Path of the clinical data sheet\n(e.g /home/user/raw_data/clinical_data.txt)")
# Control group
parser.add_argument('-ctrl', '--control', dest='control', default="NonCancer",
                	help="In the two-group comparison, which\nlabel should be considered as CONTROL\n(e.g. NonCancer)")
# Condition group
parser.add_argument('-cond', '--condition', dest='condition', default="Cancer",
                	help="In the two-group comparison, which\nlabel should be considered as CONDITION\n(e.g. Cancer)")
# Project type
parser.add_argument('-pj', '--project', dest='project', default="MyProject",
                	help="Indicate your project with one word. For example,\nin a cancer diagnosis project, you can use the\ncancer type e.g. NSCLC")
# Reference genome
parser.add_argument('-g', '--refgenome', dest='refgenome', default=referencegenome,
                	help="Path of the reference genome in fasta format\n(preferably from GENCODE)")
# Reference annotation
parser.add_argument('-ra', '--refannot', dest='refannot', default=referenceannot,
                	help="Path of the reference annotation in gtf format\n(also preferably from GENCODE)")
# Auto-detect Adapter sequence
parser.add_argument('-ad', '--autodetectadapter', dest='autodetectadapter', action="store_false",
                	help="Activate auto-detection of adapter sequences (default %(default)s)")
# Adapter sequence
parser.add_argument('-a', '--adapter', dest='adapter', metavar='', 
                	help="Path to adapter sequence (in fasta format) to be trimmed\nfrom the input data")
# Activate PhiX decontamination
parser.add_argument('-pd', '--phixdecontam', dest='phixdecontam', action='store_true',
                	help="If activated, PhiX contamination will be removed from the\ninput data (default %(default)s)")
# PhiX sequence
parser.add_argument('-p', '--phix', dest='phix', default=phix_sequences, metavar='', 
                	help="PhiX sequence to be removed from the input data")
# Allow "low" confidence in gene fusions
parser.add_argument('-l', '--low', dest='low', action="store_true", 
                	help="During gene fusion matrix construction,\nallow \"low\" confidence fusions (default %(default)s)")
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(70), metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()



# HCC 
args.project_type = "HCC"


# Main folder hosting the analysis
analysis_dir = os.path.join(parent_dir, f"{args.project_type.lower()}_analysis") if args.raw_data.endswith("data") else os.path.join(parent_dir, "validationset_analysis")

# Reporting directories
reports_dir 			  = os.path.join(analysis_dir, "reports")
pipeline_reports 		  = os.path.join(analysis_dir, "reports", "pipeline_reports")
preprocessing_reports_dir = os.path.join(analysis_dir, "reports", "preprocessing_reports")
alignment_reports_dir 	  = os.path.join(analysis_dir, "reports", "alignment_reports")
genefusion_reports_dir 	  = os.path.join(analysis_dir, "reports", "gene_fusion_reports")
varcall_reports_dir 	  = os.path.join(analysis_dir, "reports", "variant_calling_reports")
# Main sub-directories hosting the individual analysis
preprocessing_dir 		  = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir 			  = os.path.join(analysis_dir, "alignments")
transcript_quant_dir 	  = os.path.join(analysis_dir, "isoform_expression")
genefusion_dir 			  = os.path.join(analysis_dir, "gene_fusion")
expression_analysis_dir   = os.path.join(analysis_dir, "gene_expression")
altexpression_dir 		  = os.path.join(analysis_dir, "FoCT")
varcall_dir 		      = os.path.join(analysis_dir, "snv")
rnaediting_dir 			  = os.path.join(analysis_dir, "rna_editing")
final_filtmat_dir 		  = os.path.join(analysis_dir, "filtered_matrices")
prefiltmat_dir 			  = os.path.join(final_filtmat_dir, "prefiltered_matrices")
# STAR indices
star_index_dir 			  = os.path.join(analysis_dir, "alignments", "star_index")
# Temp dir of variant calling process
temp 					  = os.path.join(varcall_dir, "temp")





class data_preprocessing:

	def __init__(self):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Step1. DATA PREPROCESSING')

		global file_format
		# Creating necessary directories
		if not os.path.exists(preprocessing_dir):         os.makedirs(preprocessing_dir)  
		if not os.path.exists(preprocessing_reports_dir): os.makedirs(preprocessing_reports_dir)
		
		
		lib_type, substring, file_format = self.identify_input_libraries()
		adapter_seq = self.detect_adapters(lib_type, substring, file_format)
		self.prerocessing_and_quality_control(lib_type, substring, file_format, adapter_seq)
		self.create_star_indexes(file_format)
		return

	def identify_input_libraries(self):
		""" Detecting what type of input libraries we have. First of all 
		to what extension the raw files have and whether the input libraries
		are SE or PE """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Data preprocessing - Identifying the type of the input libraries...')

		all_file_formats = []
		files = glob.glob(os.path.join(args.raw_data, "*"))
		for _ in files:
			all_file_formats.append(_.split(".")[-1])

		extention = Counter(all_file_formats).most_common()[0][0]
		# Obtaining the raw data now that we know the extension
		data = glob.glob(os.path.join(args.raw_data, f"*.{extention}"))  ### Obtaining the raw data
		

		file_format = []
		format_extentions = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
		# FInding the exact file format of the raw data
		for random_file_names in random.sample(data, 5):
			for i, ext in enumerate(format_extentions):
				if random_file_names.endswith(ext):
					file_format.append(format_extentions[i])
		file_format = Counter(file_format).most_common()[0][0]
		
		# Checking the library type based on the names of the files
		check = []
		global lib_type
		lib_type = None
		for file_names in data:
			if file_names.endswith(file_format):
				check.append(os.path.basename(file_names).split(".")[0])
		final_counting = []
		for sample_name, occur in dict(Counter(check)).items():
			final_counting.append(occur)
		if Counter(final_counting).most_common()[0][0] == 2:
			lib_type = "PE"
		else:
			lib_type = "SE"

		substring = None
		library_substrings = [".R1.", "_R1." ,"_R1_", ".r1.", "_r1.","_r1_", "_1.", f"_1{file_format}"]
		if lib_type == "PE":
			for file_names in data:
				for j, substr in enumerate(library_substrings):
					if substr in os.path.basename(file_names):
						substring = substr
						break
		
		if lib_type == "PE":
			print(f'\t\t  ATTENTION: Auto-detection showed that the input libraries are: {lib_type} (with the sub-string: {substring}) and file extension: {file_format}')
		else:
			print(f'\t\t  ATTENTION: Auto-detection showed that the input libraries are: {lib_type} and file extension: {file_format}')
		print(f'\t\t  ATTENTION: In total {len(data)} samples will be analysed!')
		return lib_type, substring, file_format[1:]

	def detect_adapters(self, lib_type, substring, file_format):
		""" Detecting whether there are known adapters in the data """
		contains_adapters = 0
		if args.autodetectadapter:
			print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Data preprocessing - Auto-detection of adapter contamination (based on FastQC) in the input libraries...')
		
			
			temp = os.path.join(preprocessing_dir, "temp")
			
			if not os.path.exists(temp): 
				os.makedirs(temp)
				### Obtaining the raw data and we randomly pic up 10 samples to run FastQC and get the adapter sequences if any
				data = glob.glob(os.path.join(args.raw_data, f"*{file_format}"))
				if len(data) >= 10:
					data = random.sample(data, 10)
				
				startTime = datetime.now()  # Recording the start time
				runFastQC = ' '.join([
				"fastqc",  # Call fastQC to quality control all processed data
				"--threads", args.threads,  # Number of threads to use
				"--extract",  # The zipped output file will be uncompressed
				"--quiet",  # Print only log warnings
				"--outdir", temp,  # Create all output files in this specified output directory
				' '.join(data),  # String containing all samples that are about to be checked
				"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc_check_adpaterseqs-report.txt")])  # Output fastQC report
				subprocess.run(runFastQC, shell=True)
				print(f'\t\t  FastQC (with {len(data)} files) finished in {datetime.now() - startTime}')
			
			contains_adapters = []
			for sum_file in glob.glob(os.path.join(temp, "*_fastqc", "summary.txt")):
				with open(sum_file) as fin:
					for line in fin:
						if "Adapter Content" in line:
							contains_adapters.append(line.split("\t")[0])
			contains_adapters = Counter(contains_adapters).most_common()
			
			if contains_adapters[0][0] == "PASS":
				print(f'\t\t  ATTENTION: NO adapter contamination was found in the input libraries!')
				contains_adapters = None
			
			else:
				multiqc_json = os.path.join(temp, "adapter_detection_report_data", "multiqc_data.json")
				print(f'\t\t  INFO: We will now try to detect the type of adapter found in your data!')
				if not os.path.exists(multiqc_json):
					premultiqc = ' '.join([
					"multiqc",  # Call fastQC to quality control all processed data
					"--quiet",  # Print only log warnings
					"--interactive",  #  Force interactive plots
					"--title", "\'Adapter detection report\'",
					"--outdir", temp,  # Create report in the FastQC reports directory
					"--filename", "adapter_detection_report",  # Name of the output report 
					temp,  # Directory where all FastQC and Cutadapt reports reside
					"2>>", os.path.join(pipeline_reports, "preprocessing_premultiqc-report.txt")])  # Output fastQC report
					subprocess.run(premultiqc, shell=True)

				# Reading the json file from the multiqc report
				json_db = {}
				with open(multiqc_json) as json_file:
					json_db = json.load(json_file)

				# Obtaining the detected adapters from all files
				detected_adapter = []
				for content in json_db["report_plot_data"]["fastqc_adapter_content_plot"]['datasets']:
					for subcontent in content:
						detected_adapter.append(subcontent['name'].strip().split(" - ")[1])

				# Counting the adapters and the majority wins!
				adapter_name, occurances = Counter(detected_adapter).most_common()[0]
				print(f'\t\t  ATTENTION: The {adapter_name.replace("_", " ").title()} was found in your input libraries and it will be removed!')
				
				# Retrieving the adapter sequence
				adapter_file = SeqIO.parse(open(adapter_sequences), 'fasta') 
				for fasta in adapter_file:
					if fasta.id == adapter_name:
						contains_adapters = str(fasta.seq)

				# # Moving the whole auto-detecting adapters directory to the preprocessing reports directory
				# shutil.move(temp, f'{os.path.join(preprocessing_reports_dir, "autodetect_adapters")}')
		else:
			if args.adapter != None and os.path.exists(args.adapter):
				contains_adapters = args.adapter
			elif args.adapter != None and not os.path.exists(args.adapter):
				print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Data preprocessing - WARNING: The input adapter file is NOT VALID. No adapters will be removed..')
				contains_adapters = None
			else:
				print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Data preprocessing - WARNING: No custom adapter file is given and autodetection is deactivated. No adapters will be removed!')
				contains_adapters = None
		return contains_adapters

	def prerocessing_and_quality_control(self, lib_type, substring, file_format, adapter):
		""" Function that calls BBDuk to remove non-biological sequences from the input data.
		It will also remove low quality reads and short sequences. FastQC will perform a 
		thorough quality control on the preprocessed reads. 'MultiQC' will also collect 
		and  summarise  all the results in one final report (summarised_report'). 
		One can get an overview of the data before the alignment starts. """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Data preprocessing - Multilevel preprocessing of the input data...')

		os.chdir(preprocessing_dir)
		data = glob.glob(f"{args.raw_data}/*{substring.strip()}*{file_format}") if lib_type == "PE" else glob.glob(os.path.join(args.raw_data, f"*{file_format}"))  ### Obtaining the raw data


		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  BBDuk - Preprocessing, quality trimming and filtering the input data...')
		### Run BBDuk with custom settings ###
		for i, file in enumerate(data, 1):
			print(f'\t\t  {i}/{len(data)}. Preprocessing {os.path.basename(file).split(".")[0]}')
			startTime = datetime.now()  # Recording the start time
			sample_name = os.path.basename(file).split(".")[0]
			stats_file = os.path.join(preprocessing_reports_dir, f'{sample_name}.stats')
			runBBDuk = ' '.join([
			"bbduk.sh",  # Call BBDuk (BBTools) to preprocess the raw data
			f"threads={args.threads}",  # Set number of threads to use
			f"in={file}",  # Input of the forward file
			"trimpolya=20",  # Trim poly-A or poly-T tails of at least 20nt on either end of reads
			"maxns=1",  # Reads with at least 1 N (after trimming) will be discarded
			"qtrim=r",  # Trim read right ends to remove bases with quality below trimq
			"trimq=20",  # Regions with average quality BELOW this will be trimmed
			"minlength=40",  # Reads shorter than this after trimming will be discarded
			"minavgquality=20",  # Reads with average quality (after trimming) below this will be discarded (20)
			"ziplevel=7",  #  Compression level 7
			f"stats={stats_file}"])  # Write statistics about which contaminants were detected
			if adapter != None and type(adapter) == str:
				runBBDuk += ' '.join([f" literal={adapter}",
									   "ktrim=r",  # Trim reads on the right to remove bases matching reference kmers	
									   "mink=8",  # Look for shorter kmers at read tips down to this length
									   "k=12",  # Setting the kmer size we want to search for
									   "hammingdistance=1"])  # Maximum Hamming distance for ref kmers
			if args.phixdecontam and adapter != None and type(adapter) != str:
				runBBDuk += ' '.join([f" ref={args.phix},{adapter}",
									   "ktrim=r",  # Trim reads on the right to remove bases matching reference kmers	
									   "mink=8",  # Look for shorter kmers at read tips down to this length
									   "k=12",  # Setting the kmer size we want to search for
									   "hammingdistance=1"])  # Maximum Hamming distance for ref kmers
			if args.phixdecontam and adapter == None:
				runBBDuk += ' '.join([f" ref={args.phix}"])
			if lib_type == "SE":
				runBBDuk += ' '.join([f" out={preprocessing_dir}/{sample_name}.trimmed.fq.gz"])
			elif lib_type == "PE":
				runBBDuk += ' '.join([
				" tpe",
				"tbo",
				f"out={preprocessing_dir}/{sample_name}.R1.trimmed.fq.gz",  # Export filtered read to file
				f'in2={file.replace(substring,substring.replace("1","2"))}', # Import read2
				f"out2={preprocessing_dir}/{sample_name}.R2.trimmed.fq.gz"])  # Export filtered read to file
			runBBDuk += ' '.join([" 2>>", os.path.join(pipeline_reports, "preprocessing_bbduk-report.txt")])  # Output report
			subprocess.run(runBBDuk, shell=True)
			print(f'\t\t  Preprocessing finished in {datetime.now() - startTime}')

		if lib_type == "PE":
			preprocessed_R1_data = glob.glob(os.path.join(preprocessing_dir, "*R1.trimmed.fq.gz"))  ### Obtaining the raw data
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  FastQC - Quality Control (forward reads) reports are being generated...')
			runFastQCR1 = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--quiet",  # Print only log warnings
			"--outdir", preprocessing_reports_dir,  # Create all output files in this specified output directory
			' '.join(preprocessed_R1_data),  # String containing all samples that are about to be checked
			"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc_R1-report.txt")])  # Output fastQC report
			subprocess.run(runFastQCR1, shell=True)
			
			preprocessed_R2_data = glob.glob(os.path.join(preprocessing_dir, "*R2.trimmed.fq.gz"))  ### Obtaining the raw data
			print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  FastQC - Quality Control (reverse reads) reports are being generated...')
			runFastQCR2 = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--quiet",  # Print only log warnings
			"--outdir", preprocessing_reports_dir,  # Create all output files in this specified output directory
			' '.join(preprocessed_R2_data),  # String containing all samples that are about to be checked
			"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc_R2-report.txt")])  # Output fastQC report
			subprocess.run(runFastQCR2, shell=True)
		else:
			preprocessed_data = glob.glob(os.path.join(preprocessing_dir, "*.trimmed.fq.gz"))  ### Obtaining the raw data
			print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  FastQC - Quality Control reports are being generated...')
			runFastQC = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--quiet",  # Print only log warnings
			"--outdir", preprocessing_reports_dir,  # Create all output files in this specified output directory
			' '.join(preprocessed_data),  # String containing all samples that are about to be checked
			"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc-report.txt")])  # Output fastQC report
			subprocess.run(runFastQC, shell=True)

		### Run MultiQC to summarise the QC reports from all samples into a summary report
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  MultiQC - Summarising the preprocessed data...')
		runMultiQC = " ".join([
		"multiqc",  # Call MultiQC
		"--quiet",  # Print only log warnings
		"--interactive",  #  Force interactive plots
		"--title", "\'Overall preprocessing summary report\'",
		"--outdir", preprocessing_reports_dir,  # Create report in the FastQC reports directory
		"--filename", "preliminary_summarised_report",  # Name of the output report 
		preprocessing_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
		"2>>", os.path.join(pipeline_reports, "preprocessing_multiqc-report.txt")])  # Output multiQC report
		subprocess.run(runMultiQC, shell=True)

		### Moving multiqc summarised report to report directory and removing fastqc.zip files
		subprocess.run(f'mv {preprocessing_reports_dir}/preliminary_summarised_report.html {reports_dir}', shell=True)
		# subprocess.run(f'rm {preprocessing_reports_dir}/*.zip', shell=True)
		return

	def create_star_indexes(self, file_format):
		""" Building the necessary indices for the STAR aligner"""
		
		if not os.path.exists(star_index_dir): 
			os.makedirs(star_index_dir)

			star_index = " ".join([
			"STAR",  # Calling STAR software
			"--runMode genomeGenerate",  # Generate genome indexes mode
			"--runThreadN", args.threads,  # Number of cores to use for this task
			"--genomeDir", star_index_dir,  # 
			"--genomeFastaFiles", args.refgenome,  # fasta file with the ref. genome sequences
			"--sjdbGTFfile", args.refannot,  # Path to the GTF file with the annotation
			"--sjdbOverhang", self.read_length(file_format),  # Get the length of the reads
			"2>>", os.path.join(pipeline_reports, "preprocessing_STAR_genomeGenerate-report.txt")])  # Input files
			subprocess.run(star_index, shell=True)
		return

	def read_length(self, file_format):
		# Obtaining the average read length in order to build the 
		# STAR indices
		temp = os.path.join(preprocessing_dir, "temp")
		if not os.path.exists(temp):
			os.makedirs(temp)
			files = random.sample(glob.glob(os.path.join(args.raw_data, f"*{file_format}")), 4)
			runFastQC = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--extract",  # The zipped output file will be uncompressed
			"--quiet",  # Print only log warnings
			"--outdir", temp,  # Create all output files in this specified output directory
			' '.join(files),  # String containing all samples that are about to be checked
			"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc_check_adpaterseqs-report.txt")])  # Output fastQC report
			subprocess.run(runFastQC, shell=True)
		
		read_len = 0
		potential_lengths = []
		for file in os.listdir(temp):
			if file.endswith("_fastqc"):
				with open(os.path.join(temp, file, 'fastqc_data.txt')) as fin:
					for line in fin:
						if line.startswith("Sequence length"):
							detected_length = line.strip().split("\t")[1]
							if "-" in detected_length:
								detected_length = detected_length.split("-")[1]
							potential_lengths.append(detected_length)
		read_len = int(Counter(potential_lengths).most_common()[0][0])
		return str(read_len - 1)

class geneNisoform_level_analysis:
	
	def __init__(self):

		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Step2. GENE, GENE FUSION AND ISOFORM QUANTIFICATION')
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  STAR - Aligning against the reference genome...')
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Arriba - Detection of gene fusions...')
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Salmon - Estimating isoform expression...')

		if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)
		if not os.path.exists(alignment_reports_dir): os.makedirs(alignment_reports_dir)
		if not os.path.exists(expression_analysis_dir): os.makedirs(expression_analysis_dir)
		if not os.path.exists(genefusion_dir): os.makedirs(genefusion_dir)
		if not os.path.exists(genefusion_reports_dir): os.makedirs(genefusion_reports_dir)
		if not os.path.exists(transcript_quant_dir): os.makedirs(transcript_quant_dir)


		already_processed_data = [os.path.basename(files).split(".")[0] for files in glob.glob(os.path.join(alignments_dir, "*.aligned.genome.bam"))]
		preprocesed_data = glob.glob(f"{preprocessing_dir}/*.R1.trimmed.fq.gz") if lib_type == "PE" else glob.glob(f"{preprocessing_dir}/*.trimmed.fq.gz")
		for current_num, sample in enumerate(preprocesed_data, 1):
			sample_name = os.path.basename(sample).split(".")[0]
			if sample_name not in already_processed_data:
				print(f'\n\t\t  {current_num}/{len(preprocesed_data)}. Processing sample {sample_name}')
				
				### Calling the different methods
				self.aligning_against_refGenome(sample, sample_name)
				self.detect_gene_fusions(sample_name)
				self.transcript_quantification(sample_name)
			
		# Final quality check of the aligned samples
		self.mapping_quality_control()  
		return

	def aligning_against_refGenome(self, sample, sample_name):
		""" Calling STAR to align the preprocessed reads against the reference genome """		
		os.chdir(alignments_dir)  # Changing working directory cause STAR output several files in the wd
		
		outf = os.path.join(pipeline_reports, "aligning_star-detailed-report.txt")
		with open(outf, "a") as outfile:
			print(f'Aligning {sample_name}', file=outfile)
			print(f"\t\t  Aligning {sample_name}")
			startTime = datetime.now()  # Recording the start time
			STAR_align = " ".join([
			"STAR",  # Calling STAR aligner
			"--runMode alignReads",  # Map reads against the reference genome
			"--genomeDir", star_index_dir,  # Load the  ref. genome
			"--runThreadN", args.threads,  # Number of threads to be used by STAR
			"--outSAMattributes All",  # A string of desired SAM attributes, in the order desired for the output SAM
			"--sjdbOverhang", self.read_length(file_format),  # Length of the donor/acceptor sequence on each side of the junctions
			"--quantMode GeneCounts TranscriptomeSAM",  # Count reads per gene and output in transcriptome coordinates
			"--quantTranscriptomeBan IndelSoftclipSingleend",  # Prohibit indels, soft clipping and single-end alignments - compatible with RSEM
			"--twopassMode Basic",  # Basic 2-pass mapping, with all 1st pass junctions inserted into the genome indexes on the fly
			"--genomeLoad", "NoSharedMemory",  # Load genome into shared and remove it after run
			"--outSAMtype", "BAM SortedByCoordinate",  # Sort the output BAM file by Coordinate
			"--limitBAMsortRAM", "20000000000",  # Maximum available RAM (bytes) for sorting BAM (20GB)
			"--outSAMunmapped Within",  # Output the unmapped reads in the SAM file
			"--readFilesCommand gunzip -c",  # Unzipping the input .gz files
			"--chimSegmentMin 12",  # Check fusion genes
			"--chimJunctionOverhangMin 10",  # Minimum overhang for a chimeric junction
			"--chimOutJunctionFormat 0",  # Formatting type for the Chimeric.out.junction file
			"--alignSJstitchMismatchNmax 5 -1 5 5",  # Maximum number of mismatches for stitching of the splice junctions
			"--chimScoreDropMax 30",  # Max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length
			"--chimScoreJunctionNonGTAG 0",  # Penalty for a non-GT/AG chimeric junction
			"--chimScoreSeparation 1",  # Minimum difference (separation) between the best chimeric score and the next one
			"--chimSegmentReadGapMax 3",  # Maximum gap in the read sequence between chimeric segments
			"--chimMultimapNmax 50",  # Maximum number of chimeric multi-alignments
			"--chimOutType WithinBAM Junctions",  # Output old SAM into separate Chimeric.out.sam file
			"--outSAMattrRGline", f"ID:{sample_name} SM:{sample_name} PL:Illumina",  # Add SAM read group line
			"--outFileNamePrefix", f"{alignments_dir}/{sample_name}.",  # Output BAM files
			"--outReadsUnmapped Fastx", f"{alignments_dir}/{sample_name}."])  # Output of unaligned reads in fastq format
			if lib_type == "PE":
				STAR_align += ' '.join([" --readFilesIn", sample, sample.replace('R1','R2')])  # Input read
			else:
				STAR_align += ' '.join([" --readFilesIn", sample])
			STAR_align += ' '.join([" 2>>", os.path.join(pipeline_reports, "aligning_star-report.txt")])  # Input files
			# subprocess.run(STAR_align, stdout=outfile, shell=True)
			subprocess.run(STAR_align, shell=True)
			print("##############################################\n\n\n", file=outfile)

			
		subprocess.run(f'mv {alignments_dir}/{sample_name}.Aligned.sortedByCoord.out.bam {alignments_dir}/{sample_name}.aligned.genome.bam', shell=True)  # Renaming mapped (gemone) read files
		subprocess.run(f'mv {alignments_dir}/{sample_name}.Aligned.toTranscriptome.out.bam {alignments_dir}/{sample_name}.aligned.transcriptome.bam', shell=True)  # Renaming mapped (gemone) read files		
		subprocess.run(f'mv {alignments_dir}/{sample_name}.ReadsPerGene.out.tab {expression_analysis_dir}/{sample_name}.gene-expression.tsv', shell=True)  # Moving and renaming the expression files		
		if lib_type == "PE":
			subprocess.run(f'mv {alignments_dir}/{sample_name}.Unmapped.out.mate1 {alignments_dir}/{sample_name}.R1.unmapped.fastq', shell=True)  # Renaming unmapped read files
			subprocess.run(f'mv {alignments_dir}/{sample_name}.Unmapped.out.mate2 {alignments_dir}/{sample_name}.R2.unmapped.fastq', shell=True)  # Renaming unmapped read files
		else:
			subprocess.run(f'mv {alignments_dir}/{sample_name}.Unmapped.out.mate1 {alignments_dir}/{sample_name}.unmapped.fastq', shell=True)  # Renaming unmapped read files
		subprocess.run(f'pigz -9 --fast --processes {args.threads} {alignments_dir}/{sample_name}*.fastq', shell=True)  # Compress the unmapped reads
		print(f'\t\t  Alignment (STAR) finished in {datetime.now() - startTime}')
		
		### Cleaning 'alignments_dir' directory from unnecessary files
		subprocess.run(f'mv {alignments_dir}/*.out {alignment_reports_dir}', shell=True)  # Remove STAR directories
		subprocess.run(f'rm -r {alignments_dir}/{sample_name}*STAR*', shell=True)  # Remove STAR directories
		return

	def detect_gene_fusions(self, sample_name):
		""" Using the Arriba gene fusion pipeline to detect gene fusions 
		from the RNA-Seq data of all the samples. """
		outf = os.path.join(pipeline_reports, "genefusion_arriba-detailed-report.txt")
		with open(outf, "a") as outfile:
			print(f'Detecting gene fusions in {sample_name}', file=outfile)
			print(f"\t\t  Detecting gene fusions in {sample_name}")
			
			startTime = datetime.now()  # Recording the start time
			subprocess.run(f'samtools index -@ {args.threads} {alignments_dir}/{sample_name}.aligned.genome.bam', shell=True)
			runArriba = " ".join([
			arriba,  # Calling Arriba pipeline
			"-S 5",  # The 'min_support' filter discards all fusions with fewer than this many supporting reads
			"-x", f"{alignments_dir}/{sample_name}.aligned.genome.bam",  # File in BAM format with main alignments as generated by STAR (Aligned.out.bam)
			"-g", args.refannot,  # GTF file with gene annotation
			"-a", args.refgenome,  # FastA file with the reference genome sequence (assembly)
			"-b", blacklistGRCh38,  # File containing blacklisted events (recurrent artifacts and transcripts observed in healthy tissue)
			"-k", knownfusionsGRCh38,  # File containing known/recurrent fusions. Boosting sensitivity for entities characterized by fusions between the same pair of genes
			"-p", proteindomainsGRCh38,  # File in GFF3 format containing coordinates of the protein domains of genes
			"-o", f"{genefusion_dir}/{sample_name}.fusion.tsv",  # Output file with fusions that have passed all filters
			"2>>", os.path.join(pipeline_reports, "genefusion_arriba-report.txt")])  # Output report
			subprocess.run(runArriba, stdout=outfile, shell=True)
			print("##############################################\n\n\n", file=outfile)

			# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Arriba Vis (v2.0.0) - Visualisation of detected gene fusions...')
			# runArribaVis = " ".join([
			# arriba_visualisation,  # Calling Arriba visualisation pipeline
			# f"--fusions={genefusion_dir}/{sample_name}.fusion.tsv",  # Input file with all fusions that have passed all filters
			# f"--alignments={alignments_dir}/{sample_name}.aligned.genome.bam",  # File in BAM format with main alignments as generated by STAR (Aligned.out.bam)
			# f"--output={genefusion_reports_dir}/{sample_name}.fusions.pdf",  # Output as .pdf
			# f"--annotation={args.refannot}",  # GTF file with gene annotation
			# f"--cytobands={cytobandsGRCh38}",  #
			# f"--proteinDomains={proteindomainsGRCh38}",  # File in GFF3 format containing coordinates of the protein domains of genes
			# "2>>", os.path.join(pipeline_reports, "genefusion_arribaVis-report.txt")])  # Output report
			# subprocess.run(runArribaVis, shell=True)
		print(f'\t\t  Gene fusions detection (arriba) finished in {datetime.now() - startTime}')
		return

	def transcript_quantification(self, sample_name):
		""" Calling RSEM after STAR to perform isoform estimation """
		startTime = datetime.now()  # Recording the start time
		outf = os.path.join(pipeline_reports, "isoformquant_salmon-detailed-report.txt")
		with open(outf, "a") as outfile:
			print(f'Estimating isoform expression in {sample_name}', file=outfile)
			print(f"\t\t  Estimating isoform expression in {sample_name}")
			run_salmon = " ".join([
			"salmon quant",  # Calling salmon quant
			"--quiet",  # Be quiet while doing quantification
			"--useVBOpt",  # Use the Variational Bayesian EM
			"--seqBias",  # Perform sequence-specific bias correction
			"--gcBias",  # Perform fragment GC bias correction
			"--gencode",  # The transcript fasta is in GENCODE format
			"--libType U",  # Definition of the strandedness of the RNA-Seq reads
			### Bootstraps are required for estimation of technical variance
			"--numBootstraps 100",  # Number of bootstrap samples to generate
			"--threads", args.threads,  # Number of threads to be used
			"--targets", ref_salmon,  # FASTA format file containing target transcripts
			"--output", f"{transcript_quant_dir}/{sample_name}",  #  Output quantification directory
			"--alignments", f"{alignments_dir}/{sample_name}.aligned.transcriptome.bam",  # Input BAM file
			"2>>", os.path.join(pipeline_reports, "isoformquant_salmon-report.txt")])  # Input files
			subprocess.run(run_salmon, shell=True)
			print("##############################################\n\n\n", file=outfile)
			
		print(f'\t\t  Isoform expression estimation (salmon) finished in {datetime.now() - startTime}')
		return

	def mapping_quality_control(self):
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING ALIGNMENT STATS')
		### EXPORTING ALIGNMENT STATS
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC and Picard - Generating alignment stats...')
		aligned_data = glob.glob(f"{alignments_dir}/*.aligned.genome.bam")
		
		for current_num, file in enumerate(aligned_data, 1):
			file_name = os.path.basename(file).split(".")[0]
			print(f'\t\t  {current_num}/{len(aligned_data)} | Processing sample {file_name}...')
			
			qc_file = os.path.join(alignment_reports_dir, f"{file_name}.bamstat.txt")
			if not os.path.exists(qc_file):
				# RSeQC - BAM stats 
				bam_stat = ' '.join([
				"bam_stat.py",
				"-i", file,  # Input BAM file
				f"> {qc_file}",  # Output file
				"2>>", os.path.join(pipeline_reports, "postalignment_bamstats-report.txt")])
				subprocess.run(bam_stat, shell=True)

				# Picard CollectMultipleMetrics
				collectMultipleMetrics = ' '.join([
				"picard CollectMultipleMetrics",  # Call picard tools CollectMultipleMetrics
				f"INPUT= {file}",  # Input BAM file
				f"OUTPUT= {alignment_reports_dir}/{file_name}.",  # Output
				f"REFERENCE_SEQUENCE= {args.refgenome}",  # Reference sequence file
				"2>>", os.path.join(pipeline_reports, "postalignment_CollectMultipleMetrics-report.txt")])
				subprocess.run(collectMultipleMetrics, shell=True)

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  multiQC - Summarising all QC reports...')
		multiQC = " ".join([
		"multiqc",  # Call MultiQC
		"--quiet",  # Print only log warnings
		"--interactive",  #  Force interactive plots
		"--outdir", alignment_reports_dir,  # Create report in the FastQC reports directory
		"--title", "\'Overall post-alignment summary report\'",
		"--filename", "post-alignment_summarised_report",  # Name of the output report 
		alignment_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
		"2>>", os.path.join(pipeline_reports, "postalignment_mutiqc-report.txt")])  # Output multiQC report
		subprocess.run(multiQC, shell=True)

		### Moving multiqc summarised report to report directorys
		subprocess.run(f'mv {alignment_reports_dir}/*summarised_report.html {reports_dir}', shell=True)
		return

	def read_length(self, file_format):
		# Obtaining the average read length in order to build the 
		# STAR indices
		temp = os.path.join(preprocessing_dir, "temp")
		if not os.path.exists(temp):
			os.makedirs(temp)
			files = random.sample(glob.glob(os.path.join(args.raw_data, f"*{file_format}")), 4)
			runFastQC = ' '.join([
			"fastqc",  # Call fastQC to quality control all processed data
			"--threads", args.threads,  # Number of threads to use
			"--extract",  # The zipped output file will be uncompressed
			"--quiet",  # Print only log warnings
			"--outdir", temp,  # Create all output files in this specified output directory
			' '.join(files),  # String containing all samples that are about to be checked
			"2>>", os.path.join(pipeline_reports, "preprocessing_fastqc_check_adpaterseqs-report.txt")])  # Output fastQC report
			subprocess.run(runFastQC, shell=True)
		
		read_len = 0
		potential_lengths = []
		for file in os.listdir(temp):
			if file.endswith("_fastqc"):
				with open(os.path.join(temp, file, 'fastqc_data.txt')) as fin:
					for line in fin:
						if line.startswith("Sequence length"):
							detected_length = line.strip().split("\t")[1]
							if "-" in detected_length:
								detected_length = detected_length.split("-")[1]
							potential_lengths.append(detected_length)
		read_len = int(Counter(potential_lengths).most_common()[0][0])
		return str(read_len - 1)

def gene_expression_matrix():
	""" Summarisation and construction of the expression matrix on gene level. For this task each 
	alignment file (BAM) is being read and all the necessary information is being extracted. In the end 
	all info will be summasised in an expression matrix. """
	print(f'\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step3. GENERATING THE EXPRESSION MATRIX (Gene Level)')

	if not os.path.exists(prefiltmat_dir): os.makedirs(prefiltmat_dir)


	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Generating the gene-level expression matrix based on the STAR output')
	os.chdir(expression_analysis_dir)  # Changing working directory 
	files = [os.path.basename(file).split(".")[0] for file in sorted(os.listdir()) if file.endswith(".gene-expression.tsv")]
	files.insert(0, "GeneID")
	files = "\t".join(files)
	subprocess.run("paste *.gene-expression.tsv | awk 'BEGIN {OFS=\"\t\"; FS=\"\t\"}; NR>4{{j=$1; for (i=2; i<=NF; i+=4){j=j FS $i} print j}}' > gene_expression_matrix.prefiltered.tab", shell=True)
	subprocess.run(f"sed -i '1i {files}' gene_expression_matrix.prefiltered.tab", shell=True)
	subprocess.run("paste *.gene-expression.tsv | awk 'BEGIN {OFS=\"\t\"; FS=\"\t\"}; NR<5{{j=$1; for (i=2; i<=NF; i+=4){j=j FS $i} print j}}' > sample_stats.tab", shell=True)
	subprocess.run(f"sed -i '1i {files}' sample_stats.tab", shell=True)
	subprocess.run(f'mv {os.path.join(expression_analysis_dir, "gene_expression_matrix.prefiltered.tab")} {prefiltmat_dir}', shell=True)

	
	# Filtering the Gene expression matrix 
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Applying filtering step and normalisation of the data...')
	gene_expr_filt = " ".join([
	"Rscript",  # Call Rscript
	os.path.join(rscripts, "expression_filtering.R"),  # Calling the expression_filtering script
	os.path.join(prefiltmat_dir, "gene_expression_matrix.prefiltered.tab"),  # Input matrix
	args.clinical_data,  # Clinical data
	final_filtmat_dir,  # Output dir
	args.threads,  # Number of cores to use
	"gene",  # type of matrix
	f'{args.control},{args.condition}',  # Contrast for the DE
	"2>>", os.path.join(pipeline_reports, "gene_expression_filtering-report.txt")])  # Directory where all reports reside
	subprocess.run(gene_expr_filt, shell=True)
	return

def isoform_expression_matrix():
	""" Summarisation and construction of the expression matrix on gene level. For this task each 
	alignment file (BAM) is being read and all the necessary information is being extracted. In the end 
	all info will be summasised in an expression matrix. """
	print(f'\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step4. GENERATING THE EXPRESSION MATRIX (Transcript Level)')


	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Generating the isoform-level expression matrix based on the STAR output')
	isoform_quant = sorted(glob.glob(os.path.join(transcript_quant_dir, "*" + os.path.sep)))
	sample_names = [paths.split("/")[-2] for paths in isoform_quant]
	run_salmon_merge = " ".join([
	"salmon quantmerge",  # Calling salmon quantmerge
	"--quants", " ".join(isoform_quant),  # List of quantification directories
	"--names", " ".join(sample_names),  # List of names to give to the samples
	"--column numreads",  # The name of the column that will be merged; options are {len, elen, tpm, numreads}
	"--output", os.path.join(prefiltmat_dir, "isoform_expression_matrix.reads.prefiltered.tab"),  #  Output quantification directory
	"2>>", os.path.join(pipeline_reports, "isoformquant_salmon_merge-report.txt")])  # Input files
	subprocess.run(run_salmon_merge, shell=True)

	
	# Filtering the Transcript expression matrix
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Applying filtering step and normalisation of the data...')
	transcript_expr_filt = " ".join([
	"Rscript",  # Call Rscript
	os.path.join(rscripts, "expression_filtering.R"),  # Calling the expression_filtering script
	os.path.join(prefiltmat_dir, "isoform_expression_matrix.reads.prefiltered.tab"),  # Input matrix
	args.clinical_data,  # Clinical data
	final_filtmat_dir,  # Output dir
	args.threads,  # Number of cores to use
	"isoform",  # type of matrix
	f'{args.control},{args.condition}',  # Contrast for the DE
	"2>>", os.path.join(pipeline_reports, "isoform_expression_filtering-report.txt")])  # Directory where all reports reside
	subprocess.run(transcript_expr_filt, shell=True)
	return

def alternative_expression_matrix():
	""" To obtain the most abundant transcripts, 20 random healthy individuals will be picked up from 
	the clinical data file. Those 20 chosen samples will be used as a guide to elect the list of the 
	most abundant transcripts of each gene. Once we have this list, it will be used to transform the 
	transcript expression matrix into an alternative isoform expression matrix, where only the selected
	most abundant transcripts will be divided by the total tpms of all transcripts originating from the
	same gene. This way we will obtain the frequency of the alternative expression. """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step5. GENERATING THE EXPRESSION MATRIX (Alternative Isoform Expression)')


	if not os.path.exists(altexpression_dir): os.makedirs(altexpression_dir)
	json_abundtranscripts = os.path.join(altexpression_dir, "abundant_transcripts_db.json")

	
	# Obtaining the gene origin of each transcript
	database = {}
	with open(args.refannot) as refin:
		for line in refin:
			if not line.startswith("#"):
				if line.strip().split("\t")[2] == "transcript":
					database[line.strip().split("transcript_id \"")[1].split(";")[0].strip("\"")] = line.strip().split("gene_id \"")[1].split(";")[0].strip("\"")

	if not os.path.exists(json_abundtranscripts):
		most_abundant_transcripts = calc_abundant_transcripts(database)
		json.dump(most_abundant_transcripts, open(json_abundtranscripts,'w'))
	else:
		most_abundant_transcripts = json.load(open(json_abundtranscripts))
	

	# Import the transcript expression matrix with pandas as data frame
	expr_mat = pd.read_csv(os.path.join(prefiltmat_dir, 'isoform_expression_matrix.reads.prefiltered.tab'), delimiter='\t')
	# Adding the gene ID column in the data frame
	expr_mat['GeneID'] = expr_mat.apply(lambda row: database[row['Name']], axis=1)
	cols = list(expr_mat.columns.values)  # Obtaining the column names
	# Moving GeneID as first column in the data frame
	cols = [cols[-1]] + cols[0:-1]
	expr_mat = expr_mat[cols]  # Applying changes
	expr_mat = expr_mat.rename(columns={'Name': 'TranscriptID'})  # Renaming Name to transcript column
	expr_mat = expr_mat.sort_values(['GeneID','TranscriptID'], ascending=True).reset_index(drop=True)  # Sorting the data frame by GeneID and IsoformID
	# Set the index to ['GeneID', 'TranscriptID']
	expr_mat.set_index(['GeneID', 'TranscriptID'], inplace=True)


	# Obtaining the fraction of the most abundant transcript per gene 
	# divided by the sum of the transcripts belonging to that gene
	results = pd.DataFrame([])  # Saving the results to a new data frame
	for group_name, df_group in expr_mat.groupby('GeneID'):
		for index, row in df_group.iterrows():
			if index[1] in most_abundant_transcripts:  # Transcripts found in the most abundant list
				# results.append(df_group.groupby('GeneID').apply(lambda x: row/x.sum()))
				# print(df_group.groupby('GeneID').apply(lambda x: row/x.sum()))
				results = results.append(df_group.groupby('GeneID').apply(lambda x: row/x.sum()))
	# print(results)
	
	# Saving the alternative splicing matrix to a tab-delimited file
	results.round(6).to_csv(os.path.join(prefiltmat_dir, 'alternative_isoform_expression_matrix.prefiltered.tab'), na_rep='NaN', sep='\t', index=True)

	# Filtering the Alternative isoform expression matrix
	print("Applying final filtes on the alt. isoform expression matrix..")
	# Filtering the Alternative isoform expression matrix
	alt_isoform_expr_filt = " ".join([
	"Rscript",  # Call Rscript
	os.path.join(rscripts, "alt_isoform_expression_filtering.R"),  # Calling the alt_isoform_expression_filtering.R script
	os.path.join(prefiltmat_dir, "alternative_isoform_expression_matrix.prefiltered.tab"),  # Input matrix
	final_filtmat_dir,  # Output dir
	"2>>", os.path.join(pipeline_reports, "alt_isoform_expression_filtering-report.txt")])  # Directory where all reports reside
	subprocess.run(alt_isoform_expr_filt, shell=True)
	return
	
def gene_fusion_matrix():
	""" Summarisation and construction of the gene fusion expression matrix 
	based on the fusion calls that were made with the Arriba software """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step6. GENERATING THE EXPRESSION MATRIX (Gene Fusions)')


	fusion_data = glob.glob(os.path.join(genefusion_dir, "*.fusion.tsv"))

	# Obtaining all the gene fusions
	fusion_list = set()  # Keeping only unique fusions, not duplicates
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Obtaining the list of fusions detected from all samples...')
	for i, file in enumerate(fusion_data, 1):
		with open(file) as finf:
			for line in finf:
				if not line.startswith("#"):
					# Sorting the detected fusions alphabetically to avoid duplications
					fusion1 = line.strip().split("\t")[20]
					fusion2 = line.strip().split("\t")[21]
					if fusion1 != "." and fusion2 != ".":
						fusions = sorted([fusion1, fusion2])
						fusions = f"{fusions[0]}_{fusions[1]}"
						conf = line.strip().split("\t")[14]
						if not args.low and conf in ['high','medium']:  # Do now allow low confident fusions
							fusion_list.add(fusions)
						elif args.low and conf in ['high','medium','low']:  # Allow low confident fusions
							fusion_list.add(fusions)
	# print(fusion_list)
	# Using the list of the obtained fusions to make a 
	# predefined dictionary (per sample) and fill if 
	# with values while reading each file.
	overall_dict = {}
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Analysing all {len(fusion_data)} samples...')
	for i, file in enumerate(fusion_data, 1):
		sample = os.path.basename(file)[:-11]
		with open(file) as fin:
			fusion_dict = {key:0 for key in sorted(fusion_list)}
			for line in fin:
				if not line.startswith("#"):
					fusion1 = line.strip().split("\t")[20]
					fusion2 = line.strip().split("\t")[21]
					if fusion1 != "." and fusion2 != ".":
						fusions = sorted([fusion1, fusion2])
						fusions = f"{fusions[0]}_{fusions[1]}"
						conf = line.strip().split("\t")[14]
						if args.low and conf in ['high','medium','low']:
							if i==1:
								print("ATTENTION: Low confident fusions are allowed in the fusion matrix construction!")
							if fusions in fusion_dict:
								fusion_dict[fusions] += 1
						elif not args.low and conf in ['high','medium']:
							fusion_dict[fusions] += 1
			overall_dict[sample] = fusion_dict.values()

	overall_dict = {k: overall_dict[k] for k in sorted(overall_dict)}  # Sort dictionary

	# Exporting the detected fusions 
	outfile = os.path.join(prefiltmat_dir, "gene_fusion_expression_matrix.prefiltered.tab")
	with open(outfile, "w") as mat_out:
		mat_out.write("fusion\t{0}\n".format("\t".join(list(overall_dict.keys()))))
		writer = csv.writer(mat_out, delimiter='\t')
		writer.writerows(zip(sorted(fusion_list), *overall_dict.values()))

	# Filtering the Gene fusion matrix
	alt_isoform_expr_filt = " ".join([
	"Rscript",  # Call Rscript
	os.path.join(rscripts, "gene_fusion_filtering.R"),  # Calling the gene_fusion_filtering.R script
	outfile,  # Input matrix
	final_filtmat_dir,  # Output dir
	"2>>", os.path.join(pipeline_reports, "gene_fusion_filtering-report.txt")])  # Directory where all reports reside
	subprocess.run(alt_isoform_expr_filt, shell=True)
	return

class mutation_profiling:

	def __init__(self):
		
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step7. RNA EDITING and SNV PROFILING')
		
		if not os.path.exists(rnaediting_dir): os.makedirs(rnaediting_dir)
		if not os.path.exists(varcall_reports_dir): os.makedirs(varcall_reports_dir)
		
		analysed_data = []
		if os.path.exists(varcall_dir): 
			analysed_data = [os.path.basename(files).split(".")[0] for files in glob.glob(os.path.join(varcall_dir, "*.bcftools.vcf.gz"))]
		else:
			os.makedirs(varcall_dir)
		

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Bcftools - Calling somatic variants')
		aligned_data = glob.glob(f"{alignments_dir}/*.aligned.genome.bam")
		startTime_mut = datetime.now()
		for current_num, sample in enumerate(aligned_data, 1):
			sample_name = os.path.basename(sample).split(".")[0]
			if sample_name not in analysed_data:
				print(f'\n\t\t  {current_num}/{len(aligned_data)} | Processing sample {sample_name}...')
				# Calling the different methods
				self.remove_duplicates(sample, sample_name)
				self.base_recalibration(sample, sample_name)
				self.variant_calling(sample, sample_name)
				self.snv_filtering(sample, sample_name)
				self.vc_stats(sample, sample_name)
			

		# Running RediTools in a parallel fashion
		print(f'{datetime.now().strftime("%d.%m %H:%M")} 5/5 | REDItools (v2.0) - Analysing all {len(aligned_data)} samples for RNA editing events...')
		rna_edit_analysis = [os.path.basename(files).split(".")[0] for files in glob.glob(os.path.join(rnaediting_dir, "*.out.table.txt"))]
		rnaediting_args = [(sample, os.path.basename(sample).split(".")[0]) for sample in aligned_data if os.path.basename(sample).split(".")[0] not in rna_edit_analysis]
		startTime_rnaedit = datetime.now()
		rnaediting_args = [(sample, os.path.basename(sample).split(".")[0]) for sample in aligned_data]
		pool = mp.Pool(processes=int(args.threads))
		pool.map(self.RNA_editing, rnaediting_args)
		pool.close()
		print(f'\tRNA editing profiling finished after {datetime.now() - startTime_rnaedit}')

		self.rnaediting_filtering()

		print(f'\tMutation profiling finished after {datetime.now() - startTime_mut}')
		return

	def remove_duplicates(self, sample, sample_name):
		""" Removing duplicates from the aligned to the genome files """
		
		if not os.path.exists(temp):os.makedirs(temp)


		print('\t\t  1/5 | Samtools markdup - Removing duplicate reads')
		# Marking duplicate reads
		rmDuplicates = " ".join([
		"samtools rmdup",  # Call samtools markdup
		"--output-fmt BAM",  # Output in BAM format
		"-S",  # treat PE reads as SE in rmdup
		sample,
		f"{temp}/{sample_name}.rmdup.bam",
		"2>>", os.path.join(pipeline_reports, "removeDuplicates_report.txt")])
		subprocess.run(rmDuplicates, shell=True)

		# Indexing
		bam_index = " ".join([
		"samtools index",  # Call samtools markdup
		"-@", args.threads,  # Number of threads to be used
		f"{temp}/{sample_name}.rmdup.bam",
		"2>>", os.path.join(pipeline_reports, "postprocessing_bamindex-report.txt")])
		subprocess.run(bam_index, shell=True)
		return

	def base_recalibration(self, sample, sample_name):
		print('\t\t  2/5 | GATK BaseRecalibrator - Base Recalibration')
		known_sites = [f for f in glob.glob("{0}/*hg38.vcf".format(gatk_known_sites))]

		runBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk BaseRecalibrator
		"--input", f"{temp}/{sample_name}.rmdup.bam",
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--reference", args.refgenome,  # Reference sequence file
		"--output", f"{temp}/{sample_name}_recal_data.table",
		"2>>", os.path.join(pipeline_reports, "postprocessing_prebaserecalibrator-report.txt")])
		subprocess.run(runBaseRecalibrator, shell=True)

		# Generation of the recalibrated reads
		runApplyBQSR = " ".join([
		"gatk ApplyBQSR",  # Call gatk ApplyBQSR
		"--input", f"{temp}/{sample_name}.rmdup.bam",  # Input file containing sequence data
		"--bqsr-recal-file", f"{temp}/{sample_name}_recal_data.table",  # File with base quality score recalibration
		"--output", f"{temp}/{sample_name}.rmdup.recal.bam",  # Output file
		"2>>", os.path.join(pipeline_reports, "postprocessing_applybqsr-report.txt")])
		subprocess.run(runApplyBQSR, shell=True)

		# Indexing
		bam_index = " ".join([
		"samtools index",  # Call samtools markdup
		"-@", args.threads,  # Number of threads to be used
		f"{temp}/{sample_name}.rmdup.recal.bam",
		"2>>", os.path.join(pipeline_reports, "postprocessing_bamindex-report.txt")])
		subprocess.run(bam_index, shell=True)

		### Removing unnecessary files
		subprocess.run(f"rm {temp}/*.rmdup.bam", shell=True)
		subprocess.run(f"rm {temp}/*.rmdup.bam.bai", shell=True)
		subprocess.run(f"rm {temp}/*data.table", shell=True)
		return

	def variant_calling(self, sample, sample_name):
		""" Calling Bcftools to perform InDel realignment on the fly and 
		perform the  variant calling. Bcftools is a variant calling
		program for SNV, MNV, indels (<50 bp), and complex variants  """

		print('\t\t  3/5 | Bcftools - Calling somatic variants')

		## 2. Run BcfTools to call variants
		outfile = os.path.join(varcall_dir,f"{sample_name}.bcftools.vcf.gz")
		bcftools = " ".join([
		"bcftools mpileup",
		"--min-MQ 15",  # Skip alignments with mapQ smaller than 15
		"--min-BQ 15",  # Skip bases with baseQ/BAQ smaller than 15
		"--redo-BAQ",  # Recalculate BAQ on the fly
		"--per-sample-mF",  # Apply -m and -F per-sample for increased sensitivity
		"--min-ireads 2",  # Minimum number of reads for indel candidates
		"--annotate FORMAT/AD,FORMAT/DP,INFO/AD",  # Optional tags to output
		"-f", args.refgenome,  # faidx indexed reference sequence file
		"-Ou",  # Output 'u' uncompressed BCF
		os.path.join(temp, f"{sample_name}.rmdup.recal.bam"),  # Input data
		"| bcftools call",  # Calling Bcftools call
		"--threads", args.threads,  # Number of threads to use
		"--ploidy GRCh38",  # Ploidy (2)
		"--output-type z",  # Output as 'z'; compressed VCF 
		"--variants-only",  # Output variant sites only
		"--multiallelic-caller",  # alternative model for multiallelic and rare-variant calling
		"--output", outfile,
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcftools_report.txt")])
		subprocess.run(bcftools, shell=True)
		subprocess.run(f"tabix -p vcf {outfile}", shell=True)
		return

	def snv_filtering(self, sample, sample_name):
		print('\t\t  4/5 | Bcftools filter - Filtering low quality variants and removing common SNVs (based on dbSNP)')
		
		# Filtering out low quality SNPS
		outfile = os.path.join(varcall_dir, f"{sample_name}.bcftools.filt.vcf.gz")
		bcftools_filter = " ".join([
		"bcftools filter",
		"--threads", args.threads,  # Number of threads to use
		"--include", "\"QUAL>=20 && DP>=2 && FMT/DP>=2\"",  # Keep SNPs with quality >= 20, depth >= 5 reads and mapping quality >=20
		"--output-type z",  # Output as 'z'; compressed VCF 
		"--output", outfile,  # Output file
		os.path.join(varcall_dir, f"{sample_name}.bcftools.vcf.gz"),  # Input vcf
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcftools-filter_report.txt")])
		subprocess.run(bcftools_filter, shell=True)
		subprocess.run(f"tabix -p vcf {outfile}", shell=True)

		# Annotating common variants from dbSNP
		outfile = os.path.join(varcall_dir, f"{sample_name}.bcftools.filtered.annotated.vcf.gz")
		bcftools_annotate = " ".join([
		"bcftools annotate",
		"--threads", args.threads,  # Number of threads to use
		"--columns", "ID",  #
		"--annotations", common_dbsnp,  # dbSNP VCF file with annotations
		"--output-type z",  # Output as 'z'; compressed VCF 
		"--output", outfile,  # Output vcf file
		os.path.join(varcall_dir, f"{sample_name}.bcftools.filt.vcf.gz"), # Input vcf
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcftools-annotate_report.txt")])
		subprocess.run(bcftools_annotate, shell=True)
		subprocess.run(f"tabix -p vcf {outfile}", shell=True)
		os.remove(os.path.join(varcall_dir,f"{sample_name}.bcftools.filt.vcf.gz"))
		os.remove(os.path.join(varcall_dir,f"{sample_name}.bcftools.filt.vcf.gz.tbi"))

		# Exclude the common SNPs annotated from the dbSNP 
		outfile = os.path.join(varcall_dir, f"{sample_name}.bcftools.filtered.vcf.gz")
		bcftools_filter_common = " ".join([
		"bcftools filter",
		"--threads", args.threads,  # Number of threads to use
		"--include", "\"ID=='.'\"",  # Keep not annotated SNPs
		"--output-type z",  # Output as 'z'; compressed VCF 
		"--output", outfile,  # Output vcf file
		os.path.join(varcall_dir, f"{sample_name}.bcftools.filtered.annotated.vcf.gz"),  # Input vcf
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcftools-filter-common_report.txt")])
		subprocess.run(bcftools_filter_common, shell=True)
		subprocess.run(f"tabix -p vcf {outfile}", shell=True)
		os.remove(os.path.join(varcall_dir,f"{sample_name}.bcftools.filtered.annotated.vcf.gz"))
		os.remove(os.path.join(varcall_dir,f"{sample_name}.bcftools.filtered.annotated.vcf.gz.tbi"))
		return

	def vc_stats(self, sample, sample_name):
		""" Generating basic statistics and plots """
		print('\t\t  5/5 | Bcftools stats - Generating basic stats')

		bcfstats = " ".join([
		"bcftools stats",  # Calling bcfstats
		"--threads", args.threads,  # Number of threads to use for decompressing
		"--fasta-ref", args.refgenome,  # Indexed reference sequence file to determine INDEL context
		f"{varcall_dir}/{sample_name}.bcftools.filtered.vcf.gz",  # VCF file
		">", f"{varcall_reports_dir}/{sample_name}.vcf.stats.txt",
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcfstats-report.txt")])
		subprocess.run(bcfstats, shell=True)
		
		vcfplot = " ".join([
		"plot-vcfstats",  # Calling plot-vcfstats
		"--prefix", f"{varcall_reports_dir}/{sample_name}",  # Output directory/samplename
		f"{varcall_reports_dir}/{sample_name}.vcf.stats.txt",
		"2>>", os.path.join(pipeline_reports, "variantcalling_vcfplot-report.txt")])
		# subprocess.run(vcfplot, shell=True)
		return
	
	def RNA_editing(self, func_args):
		""" Generating basic statistics and plots """
		sample = func_args[0]  # Obtaining the sample
		sample_name = func_args[1]  # Obtaining the sample name
		print('\t\t  REDItools (v2.0) - Analysing the sample for RNA editing events')
		

		runReditools = " ".join([
		reditools,  # Calling REDItools
		"--file", f"{temp}/{sample_name}.rmdup.recal.bam",  # The bam file to be analyzed
		"--strict",  # Activate strict mode: only sites with edits will be included in the output
		"--strand 0",  # Strand specific protocol: 0 an unstranded protocol was followed
		"--reference", args.refgenome,  # The reference FASTA file
		"--min-edits 2",  # Minimum edit per position
		"--min-read-quality 18",  # Reads whose mapping quality is below this value will be discarded
		"--min-base-quality 15",  # Bases whose quality is below 15 will not be included in the analysis
		"--output-file", f"{rnaediting_dir}/{sample_name}.out.table.txt",  # The output statistics file
		"2>>", os.path.join(pipeline_reports, "rnaediting_reditools-report.txt")])
		subprocess.run(runReditools, shell=True)
		return 

	def rnaediting_filtering(self):
		print('\t\t  RNA editing filtering - Applying basic filters to the detected RNA editing events')

		common_snps = {}
		# Obtaining the list with the common SNPs from dbSNP database
		with open(common_dbsnp_tab) as fin:
			for line in fin:
				if not line.startswith("CHROM"):
					chrom, pos, idt, ref, alt = line.strip().split("\t")
					common_snps[f'{chrom}_{pos}_{ref}'] = alt


		rnaediting_data = glob.glob(os.path.join(rnaediting_dir,"*.out.table.txt"))
		# Performing several basic filtering steps in the detected RNA editing events		
		for i, sample in enumerate(rnaediting_data, 1):
			sample_name = os.path.basename(sample).split(".")[0]
			outfile = os.path.join(rnaediting_dir,f'{sample_name}.out.filtered.table.txt')
			with open(sample) as sin, open(outfile, 'w') as sout:
				for line in sin:
					if line.startswith("Region"):
						sout.write(line)
					else:
						coverage_q30 = int(line.strip().split("\t")[4])
						meanq = float(line.strip().split("\t")[5])
						frequency = float(line.strip().split("\t")[8])
						
						if coverage_q30 >= 10 and meanq >= 20.00 and frequency >= 0.3:
							chromosome = line.strip().split("\t")[0]
							position = line.strip().split("\t")[1]
							reference = line.strip().split("\t")[2]
							alternative = line.strip().split("\t")[7]
							common_snps_key = f'{chromosome}_{position}_{reference}'
							if len(alternative.split()) == 1:
								alternative = alternative[-1]
							else:
								alternative = alternative[0][-1]
							if common_snps_key in common_snps and alternative == common_snps[common_snps_key]:
								# print(j, common_snps_key)
								continue
							else:
								sout.write(line)
		return

def RNA_editing_matrix():
	""" Summarisation and construction of the RNA editing expression matrix 
	based on the calls that were made with the RediTools software """ 
	startTime_rnaedit = datetime.now()
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step8. GENERATING THE RNA EDITING EXPRESSION MATRIX')


	rnaediting_data = sorted(glob.glob(os.path.join(rnaediting_dir, "*.out.filtered.table.txt")))
	
	rna_edit_expression = []
	for i, file in enumerate(rnaediting_data, 1):
		samplename = os.path.basename(file)[:-23]
		# Reading each RNA editing table file
		df = pd.read_csv(file, sep='\t', usecols=['Region','Position','AllSubs'])
		df['Substitution'] = df['AllSubs'].str.split(' ').str[0]
		df = df.drop('AllSubs', axis=1)
		df = df.rename(columns={'Region': 'Chromosome'})  # Renaming 2 columns
		df = df.set_index(['Chromosome','Position','Substitution'])  # Setting the index columns
		df[samplename] = 1  # Setting 1 to all the RNA editing occurrences found in the file
		rna_edit_expression.append(df)  # Appending the recorded events to a list

	# Concatenating all per file data to a single matrix
	rna_edit_expression = pd.concat(rna_edit_expression, axis=1).fillna(0)

	# Saving the RNA-editing matrix to a file
	rna_edit_expression.to_csv(os.path.join(prefiltmat_dir, "RNAediting_expression_matrix.prefiltered.tab"), sep='\t', index=True)
	rna_edit_expression = rna_edit_expression.reset_index()  # Resetting indexes to be able to export those columns in a tsv
	rna_edit_expression.to_csv(os.path.join(varcall_dir, "RNAediting_events_to_exclude.tsv"), sep='\t', columns = ['Chromosome','Position','Substitution'], index=False)

	# Filtering out rna editing events expressed in less than 20% of the samples
	with open(os.path.join(prefiltmat_dir, "RNAediting_expression_matrix.prefiltered.tab")) as ffin, open(os.path.join(final_filtmat_dir, "RNAediting_expression_matrix.filtered.tsv"), 'w') as ffout:
		for line in ffin:
			if line.startswith("Chrom"):
				ffout.write(line)
			else:
				snp = line.strip().split("\t")
				for i, elm in enumerate(snp):
					if '0.0' in elm:
						snp[i] = '0'
					elif '1.0' in elm:
						snp[i] = '1'
				if snp[3:].count('1') >= math.ceil(((len(snp)-3)*20)/100):
					ffout.write("{0}\n".format('\t'.join(snp)))
	
	print(f'\tRNA editing expression matrix finished after {datetime.now() - startTime_rnaedit}')
	return

def snp_matrix():
	""" Summarising all SNPs found during the mutation profiling in 
	a single overall matrix (expression matrix) """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")}  Step9. GENERATING THE SNV EXPRESSION MATRIX')
	
	snp_data = ' '.join(sorted(glob.glob(f"{varcall_dir}/*.bcftools.filtered.vcf.gz")))
	

	# Merging all vcfs to one file 
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  BcfTools merge - Merging all VCF files into one and indexing...')
	subprocess.run(f'bcftools merge --output-type z --threads {args.threads} --output {os.path.join(varcall_dir,"variant_table.vcf.gz")} {snp_data}',shell = True)
	subprocess.run(f'tabix -p vcf {os.path.join(varcall_dir, "variant_table.vcf.gz")}',shell = True)  # Indexing the merged file

	# Convert the merged vcf to tab-delimited matrix
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  gatk VariantsToTable - Converting the merged VCF file to a tab-delimited file...')
	subprocess.run(f'gatk VariantsToTable -V {os.path.join(varcall_dir, "variant_table.vcf.gz")} -F CHROM -F POS -F REF -F ALT -GF GT -O {prefiltmat_dir}/snp_expression_matrix.prefiltered.tab', shell = True)  # Obtaining the GT: Genotype
	subprocess.run(f'rm {os.path.join(varcall_dir,"variant_table.vcf.*")}', shell=True)  # Removing unnecessary files
	
	snps_to_remove = []
	# Reading the RNAediting_events_to_exclude and save them to a list
	with open(os.path.join(varcall_dir, "RNAediting_events_to_exclude.tsv")) as snpin:
		for line in snpin:
			if not line.startswith("Chromosome"):
				chrom = line.strip().split("\t")[0]
				pos = line.strip().split("\t")[1]
				subs = line.strip().split("\t")[2]
				elm = f'{chrom}_{pos}_{subs}'
				snps_to_remove.append(elm)


	# Reading the SNP expression matrix and rewriting it by excluding 
	# the rows that were found in the RNA editing matrix
	# # Re-editing the SNP matrix by transforming the genotype based on the following:
	# #### SNP was not detected: 0
	# #### SNP heterozygous: 0.5
	# #### SNP homozygous: 1
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Filtering and preprocessing - Converting the genotypes in float values and performing a mild filtering step...')
	with open(os.path.join(prefiltmat_dir, "snp_expression_matrix.prefiltered.tab"), 'r') as matin, open(os.path.join(final_filtmat_dir, "snp_expression_matrix.genotype.filtered.tsv"), 'w') as matout:
		for line in matin:
			if line.startswith("CHROM"):
				matout.write(line.replace(".GT","").replace("CHROM","Chromosome").replace("POS","Position").replace("REF","Ref").replace("ALT","Alt"))
			else:
				chrmosome = line.strip().split("\t")[0]
				position = line.strip().split("\t")[1]
				ref = line.strip().split("\t")[2]
				alt = line.strip().split("\t")[3]
				element = f'{chrmosome}_{position}_{ref}{alt}'
				if not element in snps_to_remove:
					snp = line.strip().split("\t")[4:]
					snp = [elms.replace("./.","0") for elms in snp]
					for i, elm in enumerate(snp):
						if "/" in elm:
							if elm.split("/")[0] == elm.split("/")[1]:
								snp[i] = '1'
							elif elm.split("/")[0] != elm.split("/")[1]:
								snp[i] = '0.5'
					if (snp[4:].count('1') + snp[4:].count('0.5')) >= math.ceil((len(snp[4:])*20)/100):
						matout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(chrmosome, position, ref, alt,"\t".join(snp)))
	return

def calc_abundant_transcripts(database):
	""" Helper function
		Calculating the most abundant transcript per gene based on 
		20 random Healthy Control samples. We will use these transcripts
		to calculate the alternative isoform expression frequency. """
	
	random_healthyIndividuals = []  # List that will host the control group individuals
	### Reading the project's clinical data to obtain the control group individuals
	clinical_mat = pd.read_csv(args.clinical_data, sep="\t", usecols=['SampleName','Classification'])
	clinical_mat = clinical_mat.loc[clinical_mat['Classification'] == args.control]
	samples = clinical_mat['SampleName'].tolist()
	

	# Obtaining all samples from isoform quantification
	for _ in samples:
		random_healthyIndividuals.append(glob.glob(os.path.join(transcript_quant_dir, f'*{_}*', 'quant.sf'))[0])


	###  Picking up 20 individuals randomly and saving their path to random_healthyIndividuals list
	if len(random_healthyIndividuals) >= 20:
		random_healthyIndividuals = random.sample(random_healthyIndividuals, 20)
	else:
		random_healthyIndividuals = random_healthyIndividuals

	# Calculating the collective tpm per transcript per sample
	abundant_transcripts = {}
	for qunat_file in random_healthyIndividuals:
		with open(qunat_file, "r") as qin:
			for line in qin:
				if not line.startswith("Name"):
					transcript_id = str(line.strip().split("\t")[0])
					gene_id = str(database[transcript_id])
					tpm = float(line.strip().split("\t")[3])

					if gene_id in abundant_transcripts:
						if [update_value(abundant_transcripts, gene_id, idx, tpm) for idx, (isoform, value) in enumerate(abundant_transcripts[gene_id]) if transcript_id == isoform]:
							pass
						else:	
							abundant_transcripts[gene_id].append([transcript_id,tpm])
					else:
						abundant_transcripts[gene_id] = [[transcript_id,tpm]]
	

	# Obtaining the most abundant transcript from each gene 
	# and generating a list with the winner transcript.
	most_abundant_transcripts = {}
	for gene, transc_n_count in abundant_transcripts.items():
			if len(transc_n_count) == 1:
				most_abundant_transcripts[transc_n_count[0][0]] = gene
			else:
				most_abundant = max(transc_n_count, key = itemgetter(1))[0]
				most_abundant_transcripts[most_abundant] = gene
	return most_abundant_transcripts

def update_value(abundant_dict, gene_id, index, tpm):
	abundant_dict[gene_id][index][1] += tpm
	return abundant_dict

def summary_n_cleanup():
	### Run MultiQC to summarise the QC reports from all samples into a summary report
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  MultiQC - Summarising the preprocessed and aligned data...')
	runMultiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--interactive",  #  Force interactive plots
	"--title", "\'Overall summary report\'",
	"--outdir", pipeline_reports,  # Create report in the FastQC reports directory
	"--filename", "summarised_report",  # Name of the output report 
	preprocessing_reports_dir, alignment_reports_dir, varcall_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(pipeline_reports, "preprocessing_multiqc-report.txt")])  # Output multiQC report
	subprocess.run(runMultiQC, shell=True)

	### Moving multiqc summarised report to report directory and removing fastqc.zip files
	subprocess.run(f'mv {pipeline_reports}/summarised_report.html {pipeline_reports}', shell=True)

	## REMOVING UNNECESSARY FILES & REPORTS
	for path, subdir, folder in os.walk(pipeline_reports):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0:
				os.remove(file)
	
	## shutil.rmtree(temp)
	return



def main():

	# Creating reports directory
	if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)
	

	### Performing preprocessing of the data
	data_preprocessing()  # Preprocessing the raw reads

	### Performing multilevel analysis including gene and isoform level quantification
	geneNisoform_level_analysis()

	gene_expression_matrix()

	isoform_expression_matrix()

	alternative_expression_matrix()

	gene_fusion_matrix()

	### Variant profiling analysis
	mutation_profiling()

	RNA_editing_matrix()

	snp_matrix()

	summary_n_cleanup()
	

	print(f'The pipeline finished after {datetime.now() - startTime}')
	
if __name__ == "__main__": main()
