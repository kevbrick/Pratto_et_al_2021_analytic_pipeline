imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repliSim_loadModules.R')

library(data.table)

dfHiRT <- fread('hiC_v_RT.tab',header=TRUE)

## Get RT v hiC correlation Per CS
x <- flipHiCByCS(dfHiRT)

dfRet        <- plyr:::join(dfHiRT,x,by='cs')

## Flip the required chromosomes
dfRet$newHiC <- dfRet$hiC*(((dfRet$R2>0)*2)-1)

## Remove data from chromosomes with unusual eigenvector patterns
dfRet$newHiC[dfRet$cs %in% paste0('chr',c(1,14,18))] <- 0

## DF for output
dfHiC <- data.frame(hiCZyg = dfRet$newHiC)

write.table(dfHiC,
            file = 'hiCZyg.OL',
            append = FALSE,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
