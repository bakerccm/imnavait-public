# Imnavait Creek Permafrost and Active Layer Soil Microbial Communities

This repo contains code, data and partial output relating to the article:

BAKER, Christopher C.M.,  Amanda J. BARKER, Thomas A. DOUGLAS, Stacey J. DOHERTY and Robyn A. BARBATO. Seasonal variation in near-surface seasonally thawed active layer and permafrost soil microbial communities. **Environmental Research Letters** (Accepted). [doi:10.1088/1748-9326/acc542](https://doi.org/10.1088/1748-9326/acc542)

Please refer to the [published article](https://doi.org/10.1088/1748-9326/acc542) for details of our study.

## Reproducing our results

This repo **includes**:

 - [code](/code) for our bioinformatics pipeline to process the raw sequence data, based on the [dada2](https://benjjneb.github.io/dada2/tutorial.html) workflow and implemented in [Snakemake](https://snakemake.readthedocs.io);
 
 - [code](/code) for our statistical analysis in R as documented in the [snakefile](/snakefile) and summarized in the Snakemake [rulegraph](/docs) (although in reality we did not use Snakemake to execute this code);
 
 - [sample metadata](/metadata);
 
 - [soil nutrient data](/data/chemistry);
 
 - [porewater chemistry data](/data/chemistry);
 
 - [thermistor data](/data/temperature); and
 
 - [climate data](/data/climate).

This repo **does not include**:

 - our raw amplicon sequence data,
 - reference databases for taxonomic assignment, or
 - software.

To reproduce our results, this repo needs to be augmented with the raw sequence data and  reference databases that we used, and software needs to be installed.

### Sequence data

Raw amplicon sequence data files should be located in [data](/data/raw) (or a different location can be specified by updating the [config file](/config/config.yaml)). The data are available from the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under [BioProject PRJNA934356](http://www.ncbi.nlm.nih.gov/bioproject/934356). See  [data/readme.md](/data/readme.md) for further information.

### Reference databases

Reference databases for taxonomy assignment should be located in [databases](/databases) (or a different location can be specified by updating the [config file](/config/config.yaml)). For information about the reference databases, see [databases/readme.md](/databases/readme.md). Shell commands to download the data are included at [databases/download.sh](/databases/download.sh). MD5 checksums are in located in [databases/md5sums.txt](/databases/md5sums.txt).

### Software

Most software can be installed by Snakemake (or manually) using the conda environment specifications located in the [envs](/envs) folder.

The main exception is SEPP, which needs to be installed according to [instructions](https://github.com/smirarab/sepp/tree/master/sepp-package). Once installed, the [sepp_env.yaml](/envs/sepp_env.yaml) file provides a suitable environment for SEPP to run.

---
Created by [Chris Baker](https://github.com/bakerccm)
