# ELLBA (Ensemble Learning for Liquid Biopsy Analysis)

![Figure S1 GeneralScheme](https://github.com/sgiannouk/ellba/assets/22454788/7be1d871-5a02-43ca-a7b1-091b3c1e53e4)

Liquid biopsy-derived RNA sequencing (lbRNA-seq) holds immense potential for non-invasive cancer diagnostics due to its simplicity and repeatability. However, challenges arise from the lack of standardized protocols and technical variability, hampering clinical integration and leading to issues of reproducibility. To address these hurdles, we introduce a comprehensive workflow with four main objectives: 
    
    1) Harnessing the diverse molecular and functional features accessible through lbRNA-seq
    2) Implementing an Ensemble Learning framework to exploit the complementary information inherent in different feature types
    3) Rigorously benchmarking intra-sample normalization methods for improved clinical relevance
    4) Evaluating performance across independent test sets to ensure robust reproducibility. 
    
Across several datasets spanning various biosources, we observe that the optimal normalization method varies, with the straightforward Counts Per Million method performing comparably to more complex cross-sample methods. Moreover, the inclusion of novel features consistently enhances prediction accuracy beyond models based solely on gene expression. Importantly, our workflow demonstrates robust performance on entirely independent datasets, highlighting its potential for clinical applications. In summary, our approach enhances prediction accuracy and offers a promising avenue for clinical implementation, addressing vital needs in non-invasive cancer diagnostics.


INSTALLATION
------
The ELLBA software can be downloaded and installed a Docker Image that can be found here.

OUTPUT
------

USAGE
------

NOTES
------

CITATION
------

LICENSE
------
Copyright @ 2023 Stavros Giannoukakos

Licensed under the Apache License, Version 2.0


