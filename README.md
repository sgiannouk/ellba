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

OUTPUT
------


NOTES
------
- ELLBA software accepts only fastq files (accepted formats: *.fastq.gz, *.fq.gz *.fastq, *.fq)
- In case for paired end files, the read number must be indicated in both mate file name 
(e.g for read 1: ".R1.", "_R1." ,"_R1_", ".r1.", "_r1.","_r1_", or "_1.")
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

