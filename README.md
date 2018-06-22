# SPRITE4 - Pipeline for fast, scalable variant detection

This repository contains source code for SPRITE4 pipeline to perform fast and scalable variant detection. 
![workflow](workflow.png)

##Files

The folder minimap2\_parallelio\_circular contains the modified source code for [Minimap2](https://github.com/lh3/minimap2) aligner. Modification includes:

- Support for multi-node parallelism
- Multi-threading using OpenMP
- AEB, AIB output files instead of SAM format

The folder strelk2 contains the modified source code for [Strelka2](https://github.com/Illumina/strelka) variant caller. Modification includes:

- Support for multi-node parallelism
- Static load balancing to improve scalability

##Installation

The package can be installed using [BIOCONDA](https://bioconda.github.io/) recipe. The installation command is:
```
conda install -c vasupsu sprite4
```
The current recipe does not install Strelka2. We are working on this and we'll update the recipe for Strelka2 and provide test scripts and run instructions soon.

##References
1. Li, H. (2018). [Minimap2: pairwise alignment for nucleotide sequences]{https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty191/4994778}. Bioinformatics. doi:10.1093/bioinformatics/bty191
2. Kim, S., Scheffler, K. et al. (2017) [Strelka2: Fast and accurate variant calling for clinical sequencing applications](https://www.biorxiv.org/content/early/2017/09/23/192872). bioRxiv doi: 10.1101/192872

For questions, please contact Vasudevan Rengasamy, vxr162@psu.edu
