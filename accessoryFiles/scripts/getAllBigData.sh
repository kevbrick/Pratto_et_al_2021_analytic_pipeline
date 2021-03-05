#!/bin/bash

cd annotation
bash getAnnotationData.sh
cd ..

cd data/oriSSDS
bash getOriSSDSData.sh
cd ../..

cd RTSeq
bash getRDorgData.sh
cd ../../

cd RTSeq/bedgraph
bash getRTBGData.sh
cd ../../

cd RTSeq/grid
bash getGridSearchData.sh
cd ../..

cd genomeFiles
bash getGenomeFiles.sh
cd ..


