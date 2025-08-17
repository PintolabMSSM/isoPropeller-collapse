# isoPropeller-collapse
Isoform mapping and annotation stages for the isoPropeller workflow


## Installation


#### Prepare the snakemake conda environment

Installation of the required external software packages is largely handled by the pipeline itself, however a conda environment named `snakemake` needs to be present in your environment. We recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). Follow the instructions below to start the miniconda installer on Linux. When asked whether the conda environment should automatically be initialized, select 'yes'. Note that Snakemake requires the channel_priority to be set to strict. The post-installation commands to apply this setting are included in the post-installation selection below.

```bash
# Start miniconda installation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

# Post-installation commands to enforce strict channel_priority (required for Snakemake)
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

After installation and configuration of conda/miniconda, the following 'conda create' command can be used to set up the required 'snakemake' environment.

```bash
conda create -c conda-forge -c bioconda -n snakemake 'snakemake>=9.8' snakemake-executor-plugin-lsf snakemake-executor-plugin-cluster-generic 'tabulate>=0.8'
```



## Running isoPropeller collapse

Once the snakemake environment has been created, all that is needed to run the pipeline is an input file listing the the paths to csv files with long read sequencing data in bam or fastq format. For fastq input, the csv file should have two columns and include a header with the names `path_to_fastq` and `sample`, with the paths to the fastq(.gz) and sample name, respectively. Multiple fastq files can be specified per sample. These will be merged together automatically in the final analysis.

```bash
path_to_fastq,sample
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep1.fastq.gz,WTC11_CapTrap_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep2.fastq.gz,WTC11_CapTrap_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep3.fastq.gz,WTC11_CapTrap_rep3
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep1.fastq.gz,WTC11_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep2.fastq.gz,WTC11_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep3.fastq.gz,WTC11_rep3
```

The structure is similar for bam inputs, except that the header will have the name `path_to_bam` and `sample`, e.g.:

```
path_to_bam,sample
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep1.bam,WTC11_CapTrap_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep2.bam,WTC11_CapTrap_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_CapTrap_rep3.bam,WTC11_CapTrap_rep3
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep1.bam,WTC11_rep1
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep2.bam,WTC11_rep2
/sc/arion/projects/pintod02a/Xiao/NIAP_benchmark/data/WTC11_rep3.bam,WTC11_rep3
```

Once the input file is ready, run start the pipeline as follows:

```
submitjob 72 -c 4 -m 20 -q premium -P acc_pintod02b \
   ~/opt/isoPropeller-collapse/run-isoPropeller-collapse \
   -i sample-input.csv
```
