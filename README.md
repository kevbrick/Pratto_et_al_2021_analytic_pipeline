## Analytic pipeline for Pratto et al. 2020

### Pipeline accessory data (required):
The git repo contains the ONLY the pipeline script. Before running the pipeline, you need to download the accessory data. You can either run the getAccessoryFiles.sh script from the git repo:
```
bash getAccessoryFiles.sh
```
or you can get the file directly
```
wget https://hpc.nih.gov/~brickkm/meioticReplicationData/PrattoEtAl_accessoryFiles.tar.gz
tar -zxvf PrattoEtAl_accessoryFiles.tar.gz
```

### Requirements:
nextflow    : 20.01.0 \
singularity : 3.5.3 

### Configuring your nextflow environment
By default, the pipeline is configured to run from a singularity container on a SLURM-based HPC system. On such a system, it should work without reconfiguration. To run it on a local system, specify the local profile in the nextflow run command:

```
-profile singularity,local
```

Using the local profile will assume that you have a writable temporary folder named "/tmp". This, and other system-specific settings can be changed in the nextflow configuration file (accessoryFiles/config/nextflow.config.nf). 

By default, the pipeline will use up to 16 cores and 32Gb RAM. The RAM requirement is pretty inflexible, but you can change the number of cores for each process should you need to. This will simply result in a longer runtime. Resource requirements are defined in the accessoryFiles/config/nextflow.config.nf configuration file. If your system does not have the required resources, the pipeline will die once it hits the process that exceeds requirements. 

### Required environment variables
If running on a local filesystem, the environment variable TMPDIR must be set.

### Complete bash command set to run pipeline: 
```
git clone https://github.com/kevbrick/prattoEtAlAnalyticPipeline.git && cd prattoEtAlAnalyticPipeline
bash getAccessoryFiles.sh
projDir=`pwd`
nextflow run $projDir/replicationPaperPipe_v5.nf -c $projDir/accessoryFiles/config/nextflow.config.nf -profile singularity --projectdir $projDir/
```

### Other issues: 
Errors are most likely if nextflow and/or singularity are not-quite correctly configured. For example, singularity autoMounts must be enabled, and if you are using a http proxy, the proxy environment variable must be passed to singularity via the nextflow.config.nf file. Please contact me if you have questions (kevin.brick@nih.gov).  

#### ------------------------------------------------------------------------------------------------------------------------
### Alternatives to singularity container (not recommended) 
#### ------------------------------------------------------------------------------------------------------------------------
The pipeline may also be run without the singularity container. This is not recommended, however, if you insist, the easiest alternative is to use modules (accessoryFiles/modules.config.nf). 

#### Requirements (if not using the containerized pipeline: NOT RECOMMENDED): 
R/3.6.0 \
bedtools/2.27.1 \
deeptools/3.0.1 \
juicer/1.5.6 \
kallisto/0.45.0 \
macs/2.1.2 \
meme/5.0.1 \
nextflow/20.01.0 \
picard/2.9.2 \
samtools/1.9 \
sratoolkit/2.9.2 \
ucsc/388 

#### R packages: 
animation \
cowplot \
data.table \
dplyr \
extrafont \
factoextra \
flowCore \
gganimate \
ggplot2 \
ggplotify \
ggpmisc \
ggpubr \
ggridges \
gplots \
grImport2 \
grid \
gridExtra \
lsr \
Metrics \
nVennR \
numform \
optparse \
pROC \
plyr \
png \
preprocessCore \
reshape2 \
tictoc \
UpSetR \
zoo

#### Perl packages:
Bio::SeqIO (BioPerl) \
File::Basename \
File::Find \
File::Which \
Getopt::Long \
List::Util \
Math::Round \
POSIX \
Statistics::Descriptive
