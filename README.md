## Analytic pipeline for Pratto et al. 2021

### Pipeline accessory data (required):
The git repo contains the ONLY the pipeline and accessory scripts. Before running the pipeline, you need to download the accessory data. You can either run the getAccessoryFiles.sh script from the git repo:
```
bash getAccessoryFiles.sh
```
or you can get the file directly
```
wget https://hpc.nih.gov/~brickkm/meioticReplicationData/PrattoEtAl_accessoryFiles.tar.gz
tar -zxvf PrattoEtAl_accessoryFiles.tar.gz
```

### Requirements:
nextflow    : 20.10.0 \
singularity : 3.6.4 

### Configuring your nextflow environment
By default, the pipeline is configured to run from a singularity container on a SLURM-based HPC system. On such a system, it should work without reconfiguration. To run it on a local system, specify the local profile in the nextflow run command:

```
-profile singularity,local
```

Using the local profile will assume that you have a writable temporary folder named "/tmp". This, and other system-specific settings can be changed in the nextflow configuration file (accessoryFiles/config/nextflow.config.nf). 

By default, the pipeline will use up to 16 cores and 64 Gb RAM. The RAM requirement is pretty inflexible, but you can change the number of cores for each process should you need to. This will simply result in a longer runtime. Resource requirements are defined in the accessoryFiles/config/nextflow.config.nf configuration file. If your system does not have the required resources, the pipeline will die once it hits the process that exceeds requirements. 

### Required environment variables
If running on a local filesystem, the environment variable TMPDIR must be set.

### Complete bash command set to run pipeline: 
```
git clone https://github.com/kevbrick/prattoEtAlAnalyticPipeline.git && cd prattoEtAlAnalyticPipeline
bash getAccessoryFiles.sh
projDir=`pwd`
nextflow run $projDir/Pratto_et_al_2021_analyticPipeline.nf -c $projDir/accessoryFiles/config/nextflow.config.nf -profile singularity --projectdir $projDir/
```

### Other issues: 
Errors are most likely if nextflow and/or singularity are not-quite correctly configured. For example, singularity autoMounts must be enabled, and if you are using a http proxy, the proxy environment variable must be passed to singularity via the nextflow.config.nf file. Please contact me if you have questions (kevbrick@gmail.com).  

#### ------------------------------------------------------------------------------------------------------------------------
### Alternatives to singularity container (not recommended) 
#### ------------------------------------------------------------------------------------------------------------------------
The pipeline may also be run without the singularity container. This is not recommended, however, if you insist, details of all the requirements are in the singularity definition file in the top level repository folder.
