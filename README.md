# ELLBA (Ensemble Learning for Liquid Biopsy Analysis)

![Figure S1 GeneralScheme](https://github.com/sgiannouk/ellba/assets/22454788/7be1d871-5a02-43ca-a7b1-091b3c1e53e4)

Liquid biopsy-derived RNA sequencing (lbRNA-seq) holds immense potential for non-invasive cancer diagnostics due to its simplicity and repeatability. However, challenges arise from the lack of standardized protocols and technical variability, hampering clinical integration and leading to issues of reproducibility. To address these hurdles, we introduce a comprehensive workflow with four main objectives: 
    
    1) Harnessing the diverse molecular and functional features accessible through lbRNA-seq
    2) Implementing an Ensemble Learning framework to exploit the complementary information inherent 
       in different feature types
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


<p align="right">
  <svg width="107" height="32" viewBox="0 0 107 32" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M102 0H31V32H102C104.761 32 107 29.7614 107 27V5C107 2.23858 104.761 0 102 0Z" fill="#0F1418"/>
<path d="M31 0H5C2.23858 0 0 2.23858 0 5V27C0 29.7614 2.23858 32 5 32H31V0Z" fill="#366E9C"/>
<path d="M27.5592 13.0251C27.1467 11.3698 26.3646 10.1216 24.6986 10.1216H22.5505V12.6608C22.5505 14.6322 20.8791 16.2928 18.972 16.2928H13.2508C11.6865 16.2928 10.3902 17.6321 10.3902 19.2017V24.655C10.3902 26.2086 11.7401 27.1192 13.2508 27.5639C15.0614 28.0942 16.8024 28.1906 18.972 27.5639C20.413 27.146 21.8326 26.305 21.8326 24.655V22.4748H16.1168V21.7462H24.6986C26.3646 21.7462 26.9807 20.5838 27.5592 18.8427C28.1592 17.0482 28.1324 15.3232 27.5592 13.0251ZM19.3309 23.9265C19.9256 23.9265 20.4077 24.414 20.4077 25.014C20.4077 25.6193 19.9256 26.1068 19.3309 26.1068C18.7417 26.1068 18.2542 25.6139 18.2542 25.014C18.2595 24.4086 18.7417 23.9265 19.3309 23.9265ZM12.9883 15.575H18.7095C20.3005 15.575 21.5701 14.2625 21.5701 12.6662V7.20742C21.5701 5.65391 20.263 4.49144 18.7095 4.22895C16.7917 3.91289 14.7079 3.92896 12.9883 4.23431C10.5669 4.66287 10.1277 5.55748 10.1277 7.21278V9.39306H15.8543V10.1216H7.97953C6.31351 10.1216 4.85642 11.1234 4.40108 13.0251C3.8761 15.2054 3.85467 16.566 4.40108 18.8427C4.80821 20.5355 5.77782 21.7462 7.44383 21.7462H9.40984V19.132C9.40984 17.241 11.0437 15.575 12.9883 15.575ZM12.6294 7.93597C12.0347 7.93597 11.5526 7.44849 11.5526 6.84851C11.558 6.24317 12.0347 5.75569 12.6294 5.75569C13.2186 5.75569 13.7061 6.24853 13.7061 6.84851C13.7061 7.44849 13.224 7.93597 12.6294 7.93597Z" fill="white"/>
<path d="M46.8923 18.3818V22H44.5266V11.4971H48.2326C50.8791 11.4971 52.2023 12.6128 52.2023 14.8442C52.2023 15.8989 51.8215 16.7534 51.0598 17.4077C50.3029 18.0571 49.2897 18.3818 48.0202 18.3818H46.8923ZM46.8923 13.3135V16.5874H47.8225C49.0822 16.5874 49.7121 16.0356 49.7121 14.9321C49.7121 13.853 49.0822 13.3135 47.8225 13.3135H46.8923ZM61.0982 14.5L58.0513 22.6006C57.3189 24.5488 56.2154 25.5229 54.7407 25.5229C54.1792 25.5229 53.7178 25.4595 53.3565 25.3325V23.4868C53.6641 23.6675 53.9986 23.7578 54.3599 23.7578C54.9556 23.7578 55.3706 23.4771 55.605 22.9155L56.0005 21.9854L52.9536 14.5H55.5171L56.916 19.063C57.0039 19.3462 57.0723 19.6807 57.1211 20.0664H57.1504C57.1944 19.7832 57.2749 19.4536 57.3921 19.0776L58.8057 14.5H61.0982ZM66.7713 21.9121C66.4295 22.0928 65.9144 22.1831 65.2259 22.1831C63.5951 22.1831 62.7796 21.3359 62.7796 19.6416V16.2065H61.5638V14.5H62.7796V12.8813L65.0868 12.2222V14.5H66.7713V16.2065H65.0868V19.2388C65.0868 20.02 65.3968 20.4106 66.0169 20.4106C66.2611 20.4106 66.5125 20.3398 66.7713 20.1982V21.9121ZM75.4621 22H73.1549V17.7373C73.1549 16.6387 72.7546 16.0894 71.9538 16.0894C71.5436 16.0894 71.2116 16.2432 70.9577 16.5508C70.7038 16.8584 70.5768 17.249 70.5768 17.7227V22H68.2624V10.8965H70.5768V15.6133H70.6061C71.1725 14.749 71.9416 14.3169 72.9132 14.3169C74.6125 14.3169 75.4621 15.3423 75.4621 17.3931V22ZM81.0547 22.1831C79.8047 22.1831 78.8208 21.834 78.103 21.1357C77.3901 20.4326 77.0337 19.4805 77.0337 18.2793C77.0337 17.0391 77.4048 16.0698 78.147 15.3716C78.8891 14.6685 79.8926 14.3169 81.1572 14.3169C82.4023 14.3169 83.3789 14.6685 84.0869 15.3716C84.7949 16.0698 85.1489 16.9951 85.1489 18.1475C85.1489 19.3926 84.7827 20.3765 84.0503 21.0991C83.3227 21.8218 82.3242 22.1831 81.0547 22.1831ZM81.1133 16.0894C80.5664 16.0894 80.1416 16.2773 79.8388 16.6533C79.5361 17.0293 79.3847 17.5615 79.3847 18.25C79.3847 19.6904 79.9658 20.4106 81.1279 20.4106C82.2363 20.4106 82.7905 19.6709 82.7905 18.1914C82.7905 16.79 82.2314 16.0894 81.1133 16.0894ZM94.074 22H91.7669V17.8325C91.7669 16.6704 91.3519 16.0894 90.5218 16.0894C90.1214 16.0894 89.7918 16.2432 89.533 16.5508C89.2742 16.8584 89.1448 17.249 89.1448 17.7227V22H86.8304V14.5H89.1448V15.6865H89.1741C89.7259 14.7734 90.5291 14.3169 91.5838 14.3169C93.2439 14.3169 94.074 15.3472 94.074 17.4077V22Z" fill="white"/>
<path d="M102 0H5C2.23858 0 0 2.23858 0 5V27C0 29.7614 2.23858 32 5 32H102C104.761 32 107 29.7614 107 27V5C107 2.23858 104.761 0 102 0Z" fill="url(#paint0_linear)"/>
<defs>
<linearGradient id="paint0_linear" x1="0" y1="0" x2="0" y2="32" gradientUnits="userSpaceOnUse">
<stop stop-color="#BBBBBB" stop-opacity="0.1"/>
<stop offset="1" stop-opacity="0.1"/>
</linearGradient>
</defs>
</svg>
</p>
