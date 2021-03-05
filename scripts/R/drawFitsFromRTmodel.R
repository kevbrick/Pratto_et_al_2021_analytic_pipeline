#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggpubr)
library(optparse)
library(plyr)
library(png)

option_list = list(
	make_option(c("-m", "--model"),   type="character", default=NULL,  help="model file",  metavar="character"),
	make_option(c("-d", "--outdir"),  type="character", default=NULL,  help="Output Dir", metavar="character"),
	make_option(c("-o", "--outname"), type="character", default=NULL,  help="Output Name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')
source('accessoryFiles/scripts/R/simulateReplicationTiming_allCS.R')
source('accessoryFiles/scripts/R/repliSim_loadModules.R')

thisDate <- gsub(" ","_",gsub(":","",gsub("-",replacement = "",Sys.time())))

print(do.call(paste, list(names(opt), " = ", opt) ))

############################################################
imgdir  <- './'

theme7point()

load(opt$model)

if (is.null(opt$outname)){
	oNm   <- gsub('^.+/','',opt$model)
	oName <- gsub('.Rdata','',oNm)
}else{
	oName <- opt$outname
}

if (is.null(opt$outDir)){
	oDir <- paste0(getwd(),"/");
}else{
	oDir <- paste0(oDir,"/");
}

#print(paste0("*** ==> ",oName))

#print("* Making Sim BG ...")
outSimRT  <- simRT_exportRTBG(myModel = modelData$mod, myName = oName, type='sim')
#print("* Making Exp BG ...")
outRealRT <- simRT_exportRTBG(myModel = modelData$mod, myName = oName, type='exp')

print("* Moving to per-CS plots ...")

for (cs in unique(modelData$mod$model$whatCS)){
#for (cs in c("chr19")){
	print(paste0("* per-CS plots : ",cs))
	## Whole genome
	i <- simRT_drawFit(myModel     = modelData$mod,
	                   timepoint   = modelData$figData$bestStat$time,
                     pcRep       = modelData$figData$bestStat$pcRep,
								  	 myName      = oName,
								  	 chrom2use   = cs,
								  	 noNormalize = FALSE,
									   plotMe      = FALSE)

  #print(paste0('Printing fig to :',paste0(oDir,oName,'.',cs,'.png')))
  ggsave(paste0(oDir,oName,'.',cs,'.png'),i$fig ,width=6,height=2)
	#print(paste0('Printing fig to :',paste0(oDir,oName,'.',cs,'.pdf')))
  ggsave(paste0(oDir,oName,'.',cs,'.pdf'),i$fig ,width=6,height=2)
}

## Whole genome
i <- simRT_drawFit(myModel     = modelData$mod,
                   timepoint   = modelData$figData$bestStat$time,
                   pcRep       = modelData$figData$bestStat$pcRep,
		  						 myName      = oName,
			  					 chrom2use   = NULL,
				  				 noNormalize = FALSE,
								   plotMe      = FALSE)

print(paste0('Printing fig to :',paste0(oDir,oName,'.WG.png')))
ggsave(paste0(oDir,oName,'.WG.png'), i$fig, width=5, height=6)
print(paste0('Printing fig to :',paste0(oDir,oName,'.WG.pdf')))
ggsave(paste0(oDir,oName,'.WG.pdf'), i$fig, width=5, height=6)
