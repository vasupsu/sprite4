# SPRITE4 - Pipeline for fast, scalable variant detection

This repository contains source code for SPRITE4 pipeline to perform fast and scalable variant detection. 
![workflow](workflow.png)

## Assumptions and limitations

- FASTQ files are uncompressed and consist of reads < 120 bp. SPRITE4 cannot process reads > 120 bp.
- SPRITE4 can support single-end or paired-end FASTQ files. FASTQ file with interleaved paired-end reads is not supported.


## Files

The folder *sprite4\_minimap2\_modified* contains the modified source code for [Minimap2](https://github.com/lh3/minimap2) aligner. Modification includes:

- Support for multi-node parallelism using MPI
- Multi-threading using OpenMP
- AEB, AIB output files instead of SAM format
- Separate output files for different regions of reference genome

The folder *sprite4\_strelka2\_modified* contains the modified source code for [Strelka2](https://github.com/Illumina/strelka) variant caller. Modification includes:

- Support for multi-node parallelism using MPI
- AEB, AIB input file format
- Static load balancing to improve scalability

Source files included for sorting AEB, AIB files (*sampa.c*), performing SNP calling on simple reference regions (*parsnip.c*), combining the VCF output generated by strelka2 and PARSNIP in case of hybrid variant calling approach (*mergeVCF.c*)

Apart from this, test scripts (*sprite4-test, sprite4-parsnip-test, sprite4-generic-test*) are also included. 

## Requirement

- [Miniconda](https://conda.io/miniconda.html) or [Anaconda](https://www.anaconda.com/download/)
- gcc >= 4.8

## Installation

Pre-built binaries are available for Linux x86\_64 and OSX x86\_64 architectures.  The installation command is:

```
conda install -c vasupsu [-p install_path] sprite4
```

After installation, executable binaries and scripts are copied to <install\_path/bin> folder.

### Building SPRITE4 from source

Alternatively, SPRITE4 can be built from source. SPRITE4 source code is available on [GitHub](https://github.com/vasupsu/sprite4). SPRITE4 can be built using [BIOCONDA](https://bioconda.github.io) recipe located [here](https://github.com/vasupsu/bioconda-recipes/tree/master/recipes/sprite4) as below

Clone SPRITE4 BIOCONDA recipe as
```
git clone https://github.com/vasupsu/bioconda-recipes
cd bioconda-recipes/recipes
```

Build SPRITE4 as
```
conda-build sprite4 --croot <SPRITE4_build_path>
```

The built package sprite4-1.0-py27\_0.tar.bz2 an be found in the location  <SPRITE4\_build\_path/build\_platform>. Install SPRITE4 using local build as
```
conda install <SPRITE4_build_path/build_platform/sprite4-1.0-py27_0.tar.bz2>
```

## Testing using toy dataset

The [SPRITE4](https://github.com/vasupsu/sprite4) repository includes a toy dataset NA12891. Sample test script sprite4-test and sprite4-parsnip-test are simple test scripts for performing serial runs using this dataset. sprite4-test executes MAP-SAMPA-VARCALL pipeline whereas sprite4-parsnip-test executes MAP-SAMPA-PARSNIP-VARCALL pipeline shown in above Figure. Both these pipelines output the same set of variants for this dataset.These are BASH scripts and don't require any arguments.

The test script *sprite4-generic-test* has configurable parameters such as number of nodes and threads used for experiments, etc. Using the user-provided values for these parameters, the script can either execute the SPRITE pipeline or generate the commands to execute each stage of the pipeline. The latter capability is provided for multi-node runs where batch scripts are required. In this case, a user can insert the generated commands into a batch script.

## References
1. Li, H. (2018). [Minimap2: pairwise alignment for nucleotide sequences]{https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty191/4994778}. Bioinformatics. doi:10.1093/bioinformatics/bty191
2. Kim, S., Scheffler, K. et al. (2017) [Strelka2: Fast and accurate variant calling for clinical sequencing applications](https://www.biorxiv.org/content/early/2017/09/23/192872). bioRxiv doi: 10.1101/192872

For questions, please contact Vasudevan Rengasamy, vas.renga@gmail.com
