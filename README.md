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
singularity : 3.5.3 \

### Complete command set to run pipeline: 
```
git clone https://github.com/kevbrick/prattoEtAlAnalyticPipeline.git && cd prattoEtAlAnalyticPipeline
bash getAccessoryFiles.sh
projDir=`pwd`
nextflow run $projDir/replicationPaperPipe_v5.nf -c $projDir/accessoryFiles/config/nextflow.config.nf -profile singularity --projectdir $projDir/
```
## ------------------------------------------------------------------------
The pipeline may also be run without the singularity container. The easiest alternative is to use modules (accessoryFiles/modules.config.nf). \

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
