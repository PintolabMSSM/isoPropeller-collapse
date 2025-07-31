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
