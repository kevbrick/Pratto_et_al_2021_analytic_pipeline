## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

imgdir <- '.'

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(data.table)
library(nVennR)
library(UpSetR)

theme7point()

### MAIN FUNCTION
lOri <- getOriginsDFMM10('mm10_OriSSDS.Origins.tab')

allOriDF <- lOri$ori

########## SET 1: ALL ORIS in ANY WT SAMPLE : REQUIRE W/C ASYMMETRY OR BE >7Kb
oriUnionDF <- parseOriginsMM10(dfOri = allOriDF,oType='union')

########## SET 2: FINAL "HIGH CONFIDENCE" ORIS; IN >1 SAMPLE : REQUIRE W/C ASYMMETRY OR BE >7Kb
oriDF      <- parseOriginsMM10(allOriDF,oType='hiconf')

########## SCATTERPLOTS
gSc        <- plotScattersFigMM10(oriUnionDF$ori,oriDF$ori)

########## PLOT FR
gWC        <- plotWatsonCrickMM10(oriDF$ori)

########## PLOT DiNT v Oris
gDi        <- plotDiNT_v_OrisMM10(oriDF$ori)

########## PLOT DiNT v Oris
gC         <- plotOriCounts(oriDF$ori,paste0('WT_Rep',1:3))
gCU        <- plotOriCounts(oriUnionDF$ori,paste0('WT_Rep',1:3))

########## GET WT4 / RNASE SAMPLE ORIS
oriWT4     <- parseWT4Oris(allOriDF)
