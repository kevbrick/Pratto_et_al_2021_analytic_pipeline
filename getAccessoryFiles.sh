#!/bin/bash

#wget https://hpc.nih.gov/~brickkm/meioticReplicationData/PrattoEtAl_accessoryFiles.tar.gz 
# Published file is hosted on Zenodo
wget https://zenodo.org/record/4590899/files/PrattoEtAl_accessoryFiles.tar.gz

tar -zxvf PrattoEtAl_accessoryFiles.tar.gz

rm PrattoEtAl_accessoryFiles.tar.gz
