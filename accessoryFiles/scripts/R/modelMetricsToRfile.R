#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-m", "--metrics"), type="character", default=NULL,  help="model metrics file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')
source('accessoryFiles/scripts/R/simulateReplicationTiming_allCS.R')
source('accessoryFiles/scripts/R/repliSim_loadModules.R')

thisDate <- gsub(" ","_",gsub(":","",gsub("-",replacement = "",Sys.time())))

print(do.call(paste, list(names(opt), " = ", opt) ))

outFile <- gsub(pattern = 'modelMetrics.txt',replacement = 'modelMetrics.Rdata',x = opt$metrics)
theme7point()

dfMetrics <- getModelStats(opt$metrics,
                           bestModelCutoff = 0.9985)

save(dfMetrics,file=outFile)
