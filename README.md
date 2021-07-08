# SARS-CoV-2 ARTIC Panel for Illumina

This repository contains resources and information related to the [SARS-CoV-2 ARTIC Panel for Illumina](https://www.lexogen.com/sars-cov-2-whole-genome-sequencing-artic-panel/). The primers in this set are based on the published, improved sequences ([version N1](https://github.com/ItokawaK/Alt_nCov2019_primers/tree/master/Primers/ver_N1)) by [Itokawa et. al 2020](https://doi.org/10.1371/journal.pone.0239403). 

See this [this github repository](https://github.com/ItokawaK/Alt_nCov2019_primers) for more information and the original sequences.

## Repository content
The repository is split into two parts:

```
├── primers
└── example_workflow
```

1. [primers](primers)
   - Contains a file with primer sequences, binding positions and amplicion ranges with respect to the SARS-Cov2 MN908947.3 genome assembly.
3. [example_workflow](example_workflow)
   - Features a showcase script and some resources on how one could perform analysis of the obtained data. 

## Data analysis - an example workflow
As mentioned above, the [example_workflow](https://github.com/Lexogen-Tools/ARTIC_SARS-CoV-2/example_workflow) folder contains an example script as well as serveral other resources that will make it easier to get you started with your data analysis. This script implements the following steps and is based among others on tools and a workflow published by [Itokawa et. al 2020](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239403).

### The workflow
The basic steps of the analysis are outlined below. There is no need for adapter trimming if you use Illumina sequencers. This is because the shortest amplicon should be 380 bp, while Illumina sequencers only allow paired-end sequencing up to 300 bp. The workflow does include a position based primer trimming step though.
If you are using our [unique dual indicies](https://www.lexogen.com/indexing/) for barcoding also check out our [demultiplexing tool](https://github.com/Lexogen-Tools/idemuxcpp).


![analysis workflow](Lexogen_SARS-CoV-2_Workflow-Data_Analysis.png)


### Content of the workflow directory 

The example_workflow directory contains the following files or folders:
```
example_workflow
├── environment.yml
├── example_workflow.sh
├── fc_annotation.saf
└── MN908947.3
```

1. [environment.yml](example_workflow/environment.yml)
   - A yaml file to install the required dependencies into a conda environment.
   - This can easily be done by running the command ```conda env create -f environment.yml```
   - See [here](https://astrobiomike.github.io/unix/conda-intro) for an excellent tutorial on what conda is.
2. [example_workflow.sh](example_workflow/example_workflow.sh)
   - The above shown example workflow is implemented in bash
   - The script contains several variables at the beginning that might need to be edited. 
   - If you want to run this script make sure you downloaded the required scripts from [this github repository](https://github.com/ItokawaK/Alt_nCov2019_primers).
3. [fc_annotation.saf](example_workflow/fc_annotation.saf)
   - A simple annotation format file that can be used to quantify the amplicons with [featureCounts](http://subread.sourceforge.net/)
   - This file has been generated from the orginal resource published by Itokawa et. al 2020. See [here](https://doi.org/10.1371/journal.pone.0239403.s003) for the original file.
4. [MN908947.3](example_workflow/MN908947.3)
   - A folder containing the [MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) genome assembly indexed for [BWA](https://github.com/lh3/bwa) 0.7


## References

* Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M (2020) Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. PLoS ONE 15(9): [e0239403](https://doi.org/10.1371/journal.pone.0239403)

* Tools by Itokawa *et. al* used by the example workflow - https://github.com/ItokawaK/Alt_nCov2019_primers

* The original ARTIC primer panel - https://github.com/artic-network/artic-ncov2019

* SARS-Cov2 genome assembly and annotation (MN908947.3)
    * [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/MN908947)
    * [Ensembl](https://covid-19.ensembl.org/index.html)

 * [Our product page](https://www.lexogen.com/sars-cov-2-whole-genome-sequencing-artic-panel/)
