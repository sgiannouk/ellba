# ELLBA (Ensemble Learning for Liquid Biopsy Analysis)

![Figure S1 GeneralScheme](https://github.com/sgiannouk/ellba/assets/22454788/690a075f-b90a-466a-adc6-2e1e97166f13)


Liquid biopsy-derived RNA sequencing (lbRNA-seq) holds immense potential for non-invasive cancer diagnostics due to its simplicity and repeatability. However, challenges arise from the lack of standardized protocols and technical variability, hampering clinical integration and leading to issues of reproducibility. To address these hurdles, we introduce a comprehensive workflow with four main objectives: 
    
    1) Harnessing the diverse molecular and functional features accessible through lbRNA-seq
    2) Implementing an Ensemble Learning framework to exploit the complementary information inherent in different feature types
    3) Rigorously benchmarking intra-sample normalization methods for improved clinical relevance
    4) Evaluating performance across independent test sets to ensure robust reproducibility. 
    
Across several datasets spanning various biosources, we observe that the optimal normalization method varies, with the straightforward Counts Per Million method performing comparably to more complex cross-sample methods. Moreover, the inclusion of novel features consistently enhances prediction accuracy beyond models based solely on gene expression. Importantly, our workflow demonstrates robust performance on entirely independent datasets, highlighting its potential for clinical applications. In summary, our approach enhances prediction accuracy and offers a promising avenue for clinical implementation, addressing vital needs in non-invasive cancer diagnostics.


INSTALLATION
------
The ELLBA software can be downloaded and installed a Docker Image that can be found here.

USAGE
------
ELLBA operates in two distinct phases. The initial step is dedicated to the analysis of lbRNA-Seq data, involving the extraction of the six primary feature types and conducting data cleaning and normalization processes:

<pre>
usage: ellba.py [positional arguments] [optional arguments]

DESCRIPTION
An end-to-end liquid biopsy-derived RNA sequencing (lbRNA-seq) data analysis software.

positional arguments:
  -i raw_data, --inputdata raw_data
                        Path of the directory where the raw fastq or
                        fastq.gz data are residing (e.g /home/user/raw_data)
  -c clinical_data, --clindata clinical_data
                        Path of the clinical data sheet
                        (e.g /home/user/raw_data/clinical_data.txt)
  -ctrl control, --control control
                        In the two-group comparison, which
                        label should be considered as CONTROL
                        (e.g. Healthy)
  -cond condition, --condition condition
                        In the two-group comparison, which
                        label should be considered as CONDITION
                        (e.g. Cancer)
  -pj project, --project project
                        Indicate your project with one word. For example,
                        in a cancer diagnosis project, you can use the
                        cancer type e.g. NSCLC
  -g refgenome, --refgenome refgenome
                        Path of the reference genome in fasta format
                        (preferably from GENCODE)
  -ra refannot, --refannot refannot
                        Path of the reference annotation in gtf format
                        (also preferably from GENCODE)

optional arguments:
  -ad, --autodetectadapter
                        Activate auto-detection of adapter sequences (default True)
  -a , --adapter        Path to adapter sequence (in fasta format) to be trimmed
                        from the input data
  -pd, --phixdecontam   If activated, PhiX contamination will be removed from the
                        input data (default False)
  -p , --phix           PhiX sequence to be removed from the input data
  -l, --low             During gene fusion matrix construction,
                        allow "low" confidence fusions (default False)
  -t , --threads        Number of threads to be used in the analysis
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit
</pre>

<br>
<br>

During the second step, the software advances to the machine learning analysis phase. This encompasses various tasks such as data pre-processing, feature selection, and ultimately, the execution of ensemble learning techniques:

<pre>
usage: classification.py [positional arguments] [optional arguments]
    
positional arguments:
  -td trainingdir, --trainingdir trainingdir
                        Path of the directory where the
                        filtered training matrices are stored
                        (e.g /dataset/data_analysis/filtered_matrices)
  -ctrl control, --control control
                        In the two-group comparison, which
                        label should be considered as CONTROL
                        (e.g. Healthy)
  -cond condition, --condition condition
                        In the two-group comparison, which
                        label should be considered as CONDITION
                        (e.g. Cancer)
  -pj project, --project project
                        Indicate your project with one word. For example,
                        in a cancer diagnosis project, you can use the
                        cancer type e.g. NSCLC
  -ra refannot, --refannot refannot
                        Path of the reference annotation in gtf format
                        (also preferably from GENCODE)

optional arguments:
  -vd validationdir, --validationdir validationdir
                        Path of the directory where the filtered
                        external validation matrices are stored
                        (e.g /dataset/validationdata_analysis/filtered_matrices)
  -ro, --removeoutliers
                        Based on the initial exploratory analysis on the data
                        remove detected outlier samples (default False)
  -c cutoff, --cutoff cutoff
                        Minimum cutoff auc value in order for the
                        feature matrix to be utilised in the majority
                        voting classifier (default 0.65)
  -cc crosscut, --crosscut crosscut
                        Crossover cutoff threshold for the ensemble classifier (default 0.5)
  -t , --threads        Number of threads to be used in the analysis (default 5)
  -v, --version         show program's version number and exit
  -h, --help            show this help message and exit
</pre>


EXAMPLE USAGE
------
To initiate the ellba.py script, please execute the following command:
```bash
python ellba.py -i /raw/data/dir -c /raw/data/dir/clinical_data.txt -ctrl Healthy -cond Cancer -pj NSCLC -g /path/to/reference_genome.fasta -ra /path/to/reference_annotation.gtf
```

To proceed with the machine_learning.py script, please execute the following command:
```bash
python classification.py -td /dataset/data_analysis/filtered_matrices -ctrl Healthy -cond Cancer -pj NSCLC -ra /path/to/reference_annotation.gtf
```


OUTPUT
------

```bash
myproject_analysis
│
├── preprocessed_data
│   ├── SRX4472768_CLD.R1.trimmed.fq.gz
│   └── ...
├── alignments
│   ├── SRX4472768_CLD.aligned.genome.bam
│   ├── SRX4472768_CLD.aligned.transcriptome.bam
│   └── ...
├── gene_expression
│   ├── sample_stats.tab
│   ├── SRX4472768_CLD.gene-expression.tsv
│   └── ...
├── gene_fusion
│   ├── SRX4472768_CLD.fusion.tsv
│   └── ...
├── isoform_expression
│   ├── SRX4472768_CLD
│   └── ..
├── alt_isoform_expression
│   └── abundant_transcripts_db.json
├── rna_editing
│   ├── SRX4472768_CLD.out.filtered.table.txt
│   ├── SRX4472768_CLD.out.table.txt
│   └── ...
├── snv
│   ├── RNAediting_events_to_exclude.tsv
│   ├── SRX4472768_CLD.bcftools.filtered.vcf.gz
│   └── ...
├── filtered_matrices
│   ├── gene_expression_matrix.filtered.scaled.tsv
│   ├── isoform_expression_matrix.filtered.scaled.tsv
│   ├── alternative_isoform_expression_matrix.filtered.tsv
│   ├── gene_fusion_expression_matrix.filtered.tsv
│   ├── RNAediting_expression_matrix.filtered.tsv
│   ├── snp_expression_matrix.genotype.filtered.tsv
│   └── clinical_data.filtered.tsv
│   ├── plots
│   │   ├── Fig1.A.gene.pca.unfiltered.condition.jpeg
│   │   ├── Fig1.A.isoform.pca.unfiltered.condition.jpeg
│   │   └── ...
│   ├── prefiltered_matrices
│   │   ├── gene_expression_matrix.prefiltered.tab
│   │   ├── isoform_expression_matrix.reads.prefiltered.tab
│   │   ├── alternative_isoform_expression_matrix.prefiltered.tab
│   │   ├── gene_fusion_expression_matrix.prefiltered.tab
│   │   ├── RNAediting_expression_matrix.prefiltered.tab
│   │   └── snp_expression_matrix.prefiltered.tab
├── classification
│   ├── metrics_overview.txt
│   ├── selected_features.tsv
│   ├── individual_features
│   │   ├── gene_expression.ROCplot.png
│   │   ├── gene_expression.ConfusionMatrixPlot.testingset.png
│   │   └── ...
│   └── voting_classification
│       ├── gene_expression.PredResults.txt
│       ├── ...
│       ├── softVoting.TestingSet.results.tsv
│       ├── softVoting.TestingSet.ROCplot.png
│       └── softVoting.TestingSet.ConfusionMatrixPlot.png
├── reports
│   ├── preliminary_summarised_report.html
│   ├── post-alignment_summarised_report.html
│   ├── pipeline_reports
│   ├── preprocessing_reports
│   ├── alignment_reports
│   ├── gene_fusion_reports
│   └── variant_calling_reports
└──
```

NOTES
------
- ELLBA software accepts only fastq files (accepted formats: *.fastq.gz, *.fq.gz *.fastq, *.fq)
- In case for paired end files, the read number must be indicated in both mate file name <br> (e.g for read 1: ".R1.", "_R1." ,"_R1_", ".r1.", "_r1.","_r1_", or "_1.")
- 



GET HELP
------
Use the [GitHub issue tracker](https://github.com/sgiannouk/ellba/issues) to get help or to report bugs.

CITATION
------
Our article will be published soon...

LICENSE
------
Copyright @ 2023 Stavros Giannoukakos

The entire software is licensed under the Apache License, Version 2.0. <br>
The terms and conditions of both licenses can be found in the [LICENSE](https://github.com/sgiannouk/ellba/blob/main/LICENSE) file.

