# Data for Imnavait analysis

## Soil nutrient and porewater chemistry data

Soil nutrient and porewater chemistry data are provided with this repo in the [chemistry](./chemistry) folder.

## Climate data

Climate data are provided with this repo in the [climate](./climate) folder. Information about these data, including sources, is provided at [climate/readme.md](./climate/readme.md).

## Amplicon sequence data

Raw amplicon sequence data files are not supplied with this repo, and should be obtained separately before conducting analysis. Demultiplexed data are available from the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under [BioProject PRJNA934356](http://www.ncbi.nlm.nih.gov/bioproject/934356).

The demultiplexed data files should be saved to the folders [data/raw/16S](./raw/16S) and [data/raw/16S](./raw/ITS) (or another location as specified in the [config file](/config/config.yaml)) prior to conducting any analysis.

In each folder, each sample should have a read 1 file named `SAMPLE-R1.fastq.gz` and a read 2 file named `SAMPLE-R2.fastq.gz`, where `SAMPLE` is the SampleID listed in the [metadata files](/metadata).

## Thermistor data

Thermistor data are provided with this repo in the [temperature](./temperature) folder.
