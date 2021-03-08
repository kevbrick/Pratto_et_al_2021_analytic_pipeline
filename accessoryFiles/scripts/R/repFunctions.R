### ReplicationPaper functions
### KB June 12 2019
#############################

scriptFolder    <- '/accessoryFiles/scripts/R/'
imgOutputFolder <- '.'

## Fix double slashes in foldernames
scriptFolder    <- gsub('//','/',scriptFolder)
imgOutputFolder <- gsub('//','/',imgOutputFolder)

## Load bioconductor (## NOT USED CURRENTLY)
#source('http://bioconductor.org/biocLite.R')

## Load libraries
library(plyr)
library(cowplot)
library(dplyr)
library(data.table)
library(extrafont)
library(factoextra)
library(flowCore)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(gplots)
library(grid)
library(gridExtra)
library(numform) ## Format numbers nicely
library(preprocessCore)
library(reshape2)
library(rsvg)
library(tictoc)
library(png)


'%ni%'      <- Negate("%in%")

#############################
### drawOriSlice
### KB June 12 2019
### Draw a slice of data from oriSSDS experiments
### ARGS:
# SEE BELOW
### OUTPUTS:
# ggplot figure
drawOriSlice <- function(slice='chr1_38Mb_figure1',
                         slicefolder=NULL,
                         winType='w1000s147',
                         featureSize=NULL,
                         smoothingWindow = 1,
                         samples2use=NULL,
                         samples2omit=NULL,
                         pattern2omit=NULL,
                         normalizeByArea=FALSE,
                         showEnrichment=TRUE,
                         filled=FALSE,
                         revUnder=FALSE,
                         noATAC=FALSE,
                         noCGI=FALSE,
                         noG4=FALSE,
                         annotScale=.25){

  library(png)
  library(data.table)
  library(plyr)

  sliceIn <- paste0(slice,'.all.dz.tab')
  oriIn   <- paste0(slice,'.ori.bed')
  geneIn  <- paste0(slice,'.genes.bed')
  g4In    <- paste0(slice,'.g4.bed')
  CGIIn   <- paste0(slice,'.CGI.bed')
  ATACIn  <- paste0(slice,'.ATACSeq.bed')

  oriSSDSall     <- fread(sliceIn,header=TRUE)
  oriSSDSall$pos <- oriSSDSall$pos/1000000

  if (!is.null(samples2use)){
    oriSSDSall <- oriSSDSall[oriSSDSall$name %in% samples2use,]
  }

  if (!is.null(samples2omit)){
    print(paste0("Omit ",samples2omit))
    oriSSDSall <- oriSSDSall[!(oriSSDSall$name %in% samples2omit),]
  }

  if (!is.null(pattern2omit)){
    oriSSDSall <- oriSSDSall[!(grep(x = oriSSDSall$name, pattern=pattern2omit)),]
  }

  oriSSDS        <- oriSSDSall[oriSSDSall$ws==winType & oriSSDSall$fr != 'T',]

  totVal <- aggregate(oriSSDS$cover,
                      by=list(name=oriSSDS$name),
                      FUN=sum)

  names(totVal) <- c('name','totSSDS')

  meanVal <- aggregate(oriSSDS$cover,
                      by=list(name=oriSSDS$name),
                      FUN=mean)

  names(meanVal) <- c('name','meanSSDS')

  if (normalizeByArea){
    ylbl <- 'Coverage (%)'
    yNum <- '1.0'

    o <- plyr:::join(oriSSDS,totVal,by='name')
    o$cover <- o$cover / o$totSSDS

    oriSSDS <- o[,1:6]

  }else{
    if (showEnrichment){
      ylbl <- bquote('Fold enrichment')
      yNum <- '1.0'

      o <- plyr:::join(oriSSDS,meanVal,by='name')
      o$cover <- o$cover / o$meanSSDS

      oriSSDS <- o[,1:6]

    }else{
      ylbl     <- 'Coverage (FPM)'
      totVal$x <- 100
      yNum <- '1'
    }
  }

  mSSDS <- oriSSDS

  mSSDS$name <- factor(mSSDS$name,
                       levels=c(paste0('WT_Rep',1:3),
                                paste0('sperm_Rep',1:5),
                                'WT_Rep4','RNAh_WT_Rep4',
                                'WT_Rep5','RNAh_WT_Rep5',
                                'ESC_Rep1','ESC_Rep2','RNAse_ESC_Rep2'))

  xMin <- min(mSSDS$pos)
  xMax <- max(mSSDS$pos)
  
  y75 = round(max(mSSDS$cover)*.75,0)

  if (revUnder){
    mSSDS$cover[mSSDS$fr == 'R'] <- mSSDS$cover[mSSDS$fr == 'R']*-1
  }
  if (filled){
    gCov <- ggplot(mSSDS,aes(x=pos,y=cover,color=fr, fill=fr)) +
      facet_grid(name~.) +
      geom_area(lwd=.3,position='identity') +
      geom_text(aes(label=paste0("   ",name)),
                x=-Inf,
                y=Inf,
                check_overlap = TRUE,
                hjust=0,
                vjust=1,
                color='black',
                size=7*5/14 ) +
      theme(strip.background = element_blank(),
            strip.text=element_blank(),
            legend.position='none',
            axis.text.x=element_blank()) +
      xlab('') + ylab(ylbl) +
      coord_cartesian(xlim=c(xMin,xMax),expand=FALSE) +
      scale_color_manual(values=c(alpha('firebrick',.8),alpha('dodgerblue1',.8))) +
      scale_fill_manual(values=c(alpha('firebrick',.5),alpha('dodgerblue1',.5)))
    gCov
  }else{
    #geom_line(lwd=.6) +
    #geom_point(size=.3) +
    gCov <- ggplot(mSSDS,aes(x=pos,y=cover,color=fr, fill=fr)) +
      facet_grid(name~.) +
      geom_step(lwd=.3) +
      geom_text(aes(label=paste0("   ",name)),
                x=-Inf,
                y=Inf,
                check_overlap = TRUE,
                hjust=0,
                vjust=1,
                color='black',
                size=7*5/14 ) +
      theme(strip.background = element_blank(),
            strip.text=element_blank(),
            legend.position='none',
            axis.text.x=element_blank()) +
      xlab('') + ylab(ylbl) +
      scale_y_continuous(breaks=c(0,y75)) +
      coord_cartesian(xlim=c(xMin,xMax),expand=FALSE) +
      scale_fill_manual(values=c(alpha('firebrick',.1),alpha('dodgerblue1',.1))) +
      scale_color_manual(values=c('firebrick','dodgerblue1'))
  }

  ## Get genes to plot
  geneDZ <- fread(geneIn)
  names(geneDZ) <- c('cs','start','end','x1','x2','strand')

  geneDZ$from                       <- geneDZ$start
  geneDZ$from[geneDZ$strand == '-'] <- geneDZ$end[geneDZ$strand == '-']

  geneDZ$len                        <- geneDZ$end-geneDZ$start
  geneDZ$len[geneDZ$strand == '-']  <- geneDZ$len[geneDZ$strand == '-']*-1

  geneDZ$dir <- 1
  geneDZ$dir[geneDZ$strand == '-']  <- -1

  geneDZ$y <- jitter(rep(1.4,dim(geneDZ)[1]))

  ## Get oris to plot
  oriDZ <- fread(oriIn)
  names(oriDZ) <- c('cs','start','end')

  oriDZ$y <- 1.3
  
  aLbl <- data.frame(y=c(1.4,1.3),
                     lbl=c('Gene','Ori'))
  
  yFeature <- 1.3
  
  ## Get g4 to plot
  f <- file.info(g4In)$size
  if (f == 0 | noG4){
    useG4 <- FALSE
  }else{
    useG4 <- TRUE
    g4DZ <- fread(g4In)
    names(g4DZ) <- c('cs','from','to','seq','type','strand')

    g4DZ$len                      <- g4DZ$to-g4DZ$from
  
    yFeature <- yFeature-0.1
    g4DZ$y <- jitter(rep(yFeature,dim(g4DZ)[1]))
    aLbl <- rbind(aLbl,data.frame(y=yFeature,lbl=c('G4')))
  }

  ## Get ATACSeq to plot
  f <- file.info(ATACIn)$size
  if (f == 0 | noATAC){
    useATAC <- FALSE
  }else{
    useATAC <- TRUE
    atacDZ <- fread(ATACIn)
    names(atacDZ) <- c('cs','from','to')

    yFeature <- yFeature-0.1
    atacDZ$y <- jitter(rep(yFeature,dim(atacDZ)[1]))
    aLbl <- rbind(aLbl,data.frame(y=yFeature,lbl=c('ATAC')))
  }

  ## Get ATACSeq to plot
  f <- file.info(CGIIn)$size
  if (f == 0 | noCGI){
    useCGI <- FALSE
  }else{
    useCGI       <- TRUE
    cgiDZ        <- fread(CGIIn)
    names(cgiDZ) <- c('cs','from','to')

    yFeature <- yFeature-0.1
    cgiDZ$y      <- jitter(rep(yFeature,dim(cgiDZ)[1]))
    aLbl <- rbind(aLbl,data.frame(y=yFeature,lbl=c('CGI')))
  }

  cs   <- oriDZ$cs[1]

  gAnnotI <- ggplot() +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=from/1000000,
                     xend=(from+len)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.6,
                 color='grey40') +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=(from+len/2)/1000000,
                     xend=(from+len/2+dir)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.8,
                 color='grey20') +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=from/1000000,
                     xend=(from+dir)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.8,
                 color='grey20') +
    geom_rect(data=oriDZ,lwd=.3,
                 aes(xmin=start/1000000,
                     xmax=end/1000000,
                     ymin=y+0.005,
                     ymax=y-0.005),
                 size=.6, color='forestgreen') +
    coord_cartesian(xlim=c(xMin,xMax),
                    ylim=c(yFeature-0.05,1.45),
                    expand=FALSE) +
    scale_y_continuous(breaks=c(1),labels=c(yNum)) +
    theme(axis.title.y = element_text(color='white'),
          axis.ticks.y = element_blank(),
          axis.line.y=element_blank(),
          legend.position='none') +
    xlab(paste0('Position on ',cs,' (Mbp)')) +
    scale_y_continuous(breaks=aLbl$y,label=aLbl$lbl)
  #axis.text.y=element_text(color='white'),
  #    geom_label(data=aLbl,x=-Inf,hjust=0,aes(y=y,label=lbl),size=7*5/14,check_overlap = TRUE,label.size=0)+

  print(paste0("Use G4 = ",useG4))

  if (useG4){
    gAnnotI <- gAnnotI +
      geom_segment(data=g4DZ,lwd=.3,
                   aes(x=from/1000000,
                       xend=(from+len)/1000000,
                       y=y,
                       yend=y,color=strand),
                   arrow=arrow(length=unit(x=.1,units = 'cm'),
                               type='open'),
                   size=.8)
  }

  if (useATAC){
    gAnnotI <- gAnnotI +
      geom_rect(data=atacDZ,
                aes(xmin=from/1000000,
                    xmax=(to)/1000000,
                    ymin=y+0.005,
                    ymax=y-0.005),
                   color='black',
                   size=0.6)
  }

  if (useCGI){
    gAnnotI<- gAnnotI +
      geom_rect(data=cgiDZ,
                aes(xmin=from/1000000,
                    xmax=(to)/1000000,
                    ymin=y+0.005,
                    ymax=y-0.005),
                   color='orange',
                   size=0.6)
  }

  noMarg <- theme(plot.margin = unit(c(0,0,0,0),units='cm'))

  coverHeight <- 1
  annotHeight <- coverHeight*annotScale

  ## Alternative to ggarrange to align complex plots
  aP   <- cowplot::align_plots(gCov + noMarg,gAnnotI + noMarg, axis="l",align='v')
  gFig <- grid.arrange(aP[[1]],aP[[2]],
                       ncol=1,nrow=2,
                       heights=c(coverHeight,annotHeight))

  return(list(gFig    = ggarrange(gFig),
              gCover  = gCov,
              gAnnot  = gAnnotI,
              oriAll  = oriSSDSall,
              oriData = oriSSDS,
              dfAnno  = list(g4   = ifelse(useG4,g4DZ,list()),
                             atac = ifelse(useATAC,atacDZ,list()),
                             cgi  = ifelse(useCGI,cgiDZ,list()))))
}

############################
### drawOriClusterSlice
### KB June 12 2019
### Draw a slice of data from oriSSDS experiments
### ARGS:
# SEE BELOW
### OUTPUTS:
# ggplot figure
drawOriClusterSlice <- function(slice='chr1_38Mb_figure1',
                                slicefolder=NULL,
                                winType='w1000s147',
                                featureSize=NULL,
                                smoothingWindow = 1,
                                samples2use=NULL,
                                samplesOrder=NULL,
                                useOrder=FALSE,
                                samples2omit=NULL,
                                pattern2omit=NULL,
                                normalizeByArea=FALSE,
                                showEnrichment=TRUE,
                                filled=TRUE,
                                revUnder=FALSE,
                                noATAC=FALSE,
                                noCGI=FALSE,
                                annotScale=.25){

  library(png)
  library(data.table)
  library(plyr)

  if (useOrder){
    samplesOrder <- samples2use
  }

  sliceIn  <- paste0(slice,'.all.dz.tab')
  oriIn    <- paste0(slice,'.ori.bed')
  geneIn   <- paste0(slice,'.genes.bed')
  CpGpeaks <- paste0(slice,'.CpGpeaks.bed')
  CpGdens  <- paste0(slice,'.CpGdens.bedgraph')
  #g4In    <- paste0(slice,'.g4.bed')
  #CGIIn   <- paste0(slice,'.CGI.bed')
  #ATACIn  <- paste0(slice,'.ATACSeq.bed')

  oriSSDSall     <- fread(sliceIn,header=TRUE)
  oriSSDSall$pos <- oriSSDSall$pos/1000000

  if (!is.null(samples2use)){
    oriSSDSall <- oriSSDSall[oriSSDSall$name %in% samples2use,]
  }

  if (!is.null(samples2omit)){
    print(paste0("Omit ",samples2omit))
    oriSSDSall <- oriSSDSall[!(oriSSDSall$name %in% samples2omit),]
  }

  if (!is.null(pattern2omit)){
    oriSSDSall <- oriSSDSall[!(grep(x = oriSSDSall$name, pattern=pattern2omit)),]
  }

  oriSSDS        <- oriSSDSall[oriSSDSall$ws==winType & oriSSDSall$fr != 'T',]

  totVal <- aggregate(oriSSDS$cover,
                      by=list(name=oriSSDS$name),
                      FUN=sum)

  names(totVal) <- c('name','totSSDS')

  meanVal <- aggregate(oriSSDS$cover,
                       by=list(name=oriSSDS$name),
                       FUN=mean)

  names(meanVal) <- c('name','meanSSDS')

  if (normalizeByArea){
    ylbl <- 'Coverage (%)'
    yNum <- '1.0'

    o <- plyr:::join(oriSSDS,totVal,by='name')
    o$cover <- o$cover / o$totSSDS

    oriSSDS <- o[,1:6]

  }else{
    if (showEnrichment){
      ylbl <- bquote('Ori-SSDS Fold enrichment')
      yNum <- '1.0'

      o <- plyr:::join(oriSSDS,meanVal,by='name')
      o$cover <- o$cover / o$meanSSDS

      oriSSDS <- o[,1:6]

    }else{
      ylbl     <- 'Coverage (FPM)'
      totVal$x <- 100
      yNum <- '1'
    }
  }

  #print(dim(oriSSDS))
  #print(head(oriSSDS))

  mSSDS <- oriSSDS

  if (!is.null(samplesOrder)){
    mSSDS$name <- factor(mSSDS$name, levels=c(samplesOrder))
  }else{

    mSSDS$name <- factor(mSSDS$name,
                         levels=c(paste0('WT_Rep',1:3),
                                  paste0('sperm_Rep',1:5),
                                  paste0('WT_Rep',4:100),
                                  paste0('RNAh_WT_Rep',1:100)))
  }

  xMin <- min(mSSDS$pos)
  xMax <- max(mSSDS$pos)

  if (revUnder){
    mSSDS$cover[mSSDS$fr == 'R'] <- mSSDS$cover[mSSDS$fr == 'R']*-1
  }

  ## Get y ticks
  nClose <- 10^floor(log10(max(mSSDS$cover)))
  yTicks  <- c(0,round(max(mSSDS$cover)/nClose*.75) * nClose)

  if (filled){
    gCov <- ggplot(mSSDS,aes(x=pos,y=cover,color=fr, fill=fr)) +
      facet_grid(name~.) +
      geom_area(lwd=.3,position='identity') +
      geom_text(aes(label=paste0("   ",name)),
                x=-Inf,
                y=Inf,
                check_overlap = TRUE,
                hjust=0,
                vjust=1,
                color='black',
                size=7*5/14 ) +
      theme(strip.background = element_blank(),
            strip.text=element_blank(),
            legend.position='none',
            axis.text.x=element_blank()) +
      xlab('') + ylab(ylbl) +
      coord_cartesian(xlim=c(xMin,xMax),expand=FALSE) +
      scale_color_manual(values=c(alpha('firebrick',.8),alpha('dodgerblue1',.8))) +
      scale_fill_manual(values=c(alpha('firebrick',.5),alpha('dodgerblue1',.5))) +
      scale_y_continuous(breaks=yTicks)
    gCov
  }else{
    #geom_line(lwd=.6) +
    #geom_point(size=.3) +
    gCov <- ggplot(mSSDS,aes(x=pos,y=cover/x*100,color=fr, fill=fr)) +
      facet_grid(name~.) +
      geom_step(lwd=.5) +
      geom_text(aes(label=paste0("   ",name)),
                x=-Inf,
                y=Inf,
                check_overlap = TRUE,
                hjust=0,
                vjust=1,
                color='black',
                size=7*5/14 ) +
      theme(strip.background = element_blank(),
            strip.text=element_blank(),
            legend.position='none',
            axis.text.x=element_blank()) +
      xlab('') + ylab(ylbl) +
      coord_cartesian(xlim=c(xMin,xMax),expand=FALSE) +
      scale_fill_manual(values=c(alpha('firebrick',.1),alpha('dodgerblue1',.1))) +
      scale_color_manual(values=c('firebrick','dodgerblue1'))
  }

  ## Get genes to plot
  geneDZ <- fread(geneIn)
  names(geneDZ) <- c('cs','start','end','x1','x2','strand')

  geneDZ$from                       <- geneDZ$start
  geneDZ$from[geneDZ$strand == '-'] <- geneDZ$end[geneDZ$strand == '-']

  geneDZ$len                        <- geneDZ$end-geneDZ$start
  geneDZ$len[geneDZ$strand == '-']  <- geneDZ$len[geneDZ$strand == '-']*-1

  geneDZ$dir <- 1
  geneDZ$dir[geneDZ$strand == '-']  <- -1

  geneDZ$y <- jitter(rep(1.3,dim(geneDZ)[1]))

  ## Get oris to plot
  oriDZ <- fread(oriIn)
  names(oriDZ) <- c('cs','start','end')

  oriDZ$y <- 1.4

  ## Get CpG peaks to plot
  f <- file.info(CpGpeaks)$size
  if (f == 0 || noCGI){
    useCpGpk     <- FALSE
  }else{
    useCpGpk     <- TRUE
    cgiDZ        <- fread(CpGpeaks)
    names(cgiDZ) <- c('cs','from','to')

    #cgiDZ$y     <- jitter(rep(1.0,dim(cgiDZ)[1]))
    cgiDZ$y      <- rep(1.0,dim(cgiDZ)[1])
  }

  aLbl <- data.frame(y=c(1.4,1.3,1.2,1.1,1.0),
                     lbl=c('Ori','Gene','G4s','ATAC','CGI'))
  #lbl=c('Transcripts','Origins','G4s','ATAC-Seq','CGI'))

  cs   <- oriDZ$cs[1]

  gAnnotI <- ggplot() +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=from/1000000,
                     xend=(from+len)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.6,
                 color='grey40') +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=(from+len/2)/1000000,
                     xend=(from+len/2+dir)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.8,
                 color='grey20') +
    geom_segment(data=geneDZ,lwd=.3,
                 aes(x=from/1000000,
                     xend=(from+dir)/1000000,
                     y=y,
                     yend=y),
                 arrow=arrow(length=unit(x=.1,units = 'cm'),
                             type='open'),
                 size=.8,
                 color='grey20') +
    geom_rect(data=oriDZ,lwd=.3,
              aes(xmin=start/1000000,
                  xmax=end/1000000,
                  ymin=y+0.005,
                  ymax=y-0.005),
              size=.6, fill='firebrick',
              color='orange') +
    coord_cartesian(xlim=c(xMin,xMax),
                    ylim=c(1.25,1.45),
                    expand=FALSE) +
    scale_y_continuous(breaks=c(1),labels=c(yNum)) +
    theme(axis.title.y = element_text(color='white'),
          axis.ticks.y = element_blank(),
          axis.line.y=element_blank(),
          legend.position='none') +
    xlab(paste0('Position on ',cs,' (Mbp); Region width = ',(round((xMax-xMin)*1000000/1000/10)*10),' Kb')) +
    scale_y_continuous(breaks=aLbl$y,label=aLbl$lbl)
  #axis.text.y=element_text(color='white'),
  #    geom_label(data=aLbl,x=-Inf,hjust=0,aes(y=y,label=lbl),size=7*5/14,check_overlap = TRUE,label.size=0)+

  ### CpG density plot
  if (useCpGpk){
    cgp <- read.table(CpGpeaks)
    names(cgp) <- c('cs','from','to')
    cgp$from   <- cgp$from / 1000000
    cgp$to     <- cgp$to   / 1000000
  }

  cgd <- read.table(CpGdens)
  names(cgd) <- c('cs','from','to','cpg')
  cgd$from   <- cgd$from / 1000000
  cgd$to     <- cgd$to   / 1000000
  cgd$cpgPC  <- cgd$cpg   / 1000*100 #win size = 1kb
  cgMax      <- max(cgd$cpgPC)

  ## Get CpG y ticks
  nCGClose  <- 10^floor(log10(cgMax))
  yCGTicks  <- c(0,round(cgMax/nCGClose*.75) * nCGClose)

  gCpG <- ggplot(cgd) +
    geom_line(aes(x=from,y=cpgPC)) +
    geom_area(aes(x=from,y=cpgPC),
              fill='grey70') +
    theme(strip.background = element_blank(),
          strip.text=element_blank(),
          legend.position='none',
          axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab('') +
    ylab('CpG (%)') +
    scale_y_continuous(breaks=c(0,yCGTicks)) +
    coord_cartesian(xlim=c(xMin,xMax),
                    ylim=c(-1*(cgMax/8),cgMax),
                    expand=FALSE)

  if (useCpGpk){
    gCpG <- gCpG + geom_rect(data=cgp,
              aes(xmin=(from-500/1000000),
                  xmax=(to  +500/1000000)),
              ymin=-1*(cgMax/10),
              ymax=-1*(cgMax/60),
              lwd=4,
              fill='grey70')

  }

  noMarg   <- theme(plot.margin = unit(c(0,0,0,0.05),units='cm'))
  negMarg  <- theme(plot.margin = unit(c(-.1,0,-0.1,0.05),units='cm'))
  bottMarg <- theme(plot.margin = unit(c(-.1,0,0.05,0.05),units='cm'))

  coverHeight <- 1
  annotHeight <- coverHeight*annotScale

  ## Alternative to ggarrange to align complex plots
  aP   <- cowplot::align_plots(gCov + noMarg,
                      gCpG + negMarg,
                      gAnnotI + bottMarg,
                      axis="l",align='v')

  gFig <- grid.arrange(aP[[1]],aP[[2]],aP[[3]],
                       ncol=1,nrow=3,
                       heights=c(coverHeight,annotHeight*1.25,annotHeight))

  return(gFig    = ggarrange(gFig))
}

#############################
### drawRTvRepSlice
### KB Oct 12 2019
### Draw a slice of DSBs v RT
### ARGS:
# SEE BELOW
### OUTPUTS:
# ggplot figure
drawRTvRepSlice <- function(slice='chr12',
                         slicefolder='.',
                         winType='w150ks50k',
                         smoothingWindow = 1,
                         normalizeByArea=FALSE,
                         filled=TRUE,
                         annotScale=.25){

  library(png)
  library(data.table)
  library(plyr)

  print(slice)
  sliceIn <- paste0(slice,'.dz.tab')

  sAll     <- fread(sliceIn,header=TRUE)
  sAll$pos <- sAll$pos/1000000

  dDSB    <- sAll[sAll$name == 'spo11',]
  dfRT     <- sAll[sAll$name == 'RT',]

  dfRT <- dfRT[seq(1,length(dfRT$pos),by=5),]

  dfDSB    <- data.frame(pos   = rollapply(dDSB$pos,7,function(x){sum(x)/length(x)}),
                         cover = rollapply(dDSB$cover,7,function(x){sum(x)/length(x)}),
                         name='spo11')

  noLeg    <- theme(legend.position='none')
  noX      <- theme(axis.text.x=element_blank())
  noMarg   <- theme(plot.margin = unit(c(0,0,0,0),units='cm'))

  maxX     <- ceiling(max(sAll$pos)/5)*5

  dfRT <- dfRT[dfRT$pos > 3 & dfRT$pos < 120.12,]

  dfDSB$cover <- standardizeMNSD(dfDSB$cover) - min(standardizeMNSD(dfDSB$cover))
  dfRT$cover <- standardizeMNSD(dfRT$cover) - min(standardizeMNSD(dfRT$cover))

#  dfX2 <- rbind(dfDSB,dfRT)
#    geom_bar(data=dfDSB,
  #aes(x=pos,y=cover,color=name),
  #stat='identity',color='grey30',lwd=.2) +
  gX2 <- ggplot() +
    geom_area(data=dfRT,
              aes(x=pos,y=cover,color=name),
              lwd=.2,fill='firebrick',alpha=.1) +
    geom_line(data=dfDSB,
              aes(x=pos,y=cover,color=name),
              stat='identity',color='grey30',lwd=.4) +
    geom_line(data=dfRT,
              aes(x=pos,y=cover,color=name),
              lwd=.2,color='firebrick') +
    xlab('Position on chr12 (Mb)') +
    ylab('Density') +
    coord_cartesian(xlim=c(0,maxX),expand=FALSE) +
    noLeg + noMarg  +
    scale_x_continuous(breaks=seq(0,maxX,by=50)) +
    scale_y_continuous(breaks=c(0,9),labels=c(0,9)) +
    annotate(geom='text',x=maxX,y=Inf,label='DSBs (Spo11)  ',color='grey30'   ,hjust=1,vjust=1,size=7*5/14) +
    annotate(geom='text',x=maxX,y=Inf,label='\nMeiotic RT  '  ,color='firebrick',hjust=1,vjust=1,size=7*5/14)
#
#   gDSB    <- ggplot(dfDSB,
#                      aes(x=pos,y=cover,color=name)) +
#     geom_bar(stat='identity',color='grey40',lwd=.2) +
#     xlab('') +
#     ylab('DSBs (Spo11)') +
#     coord_cartesian(xlim=c(0,maxX),expand=FALSE) +
#     noLeg + noX + noMarg
#
#   gRT       <- ggplot(dfRT,
#                       aes(x=pos,y=cover*-1,color=name)) +
#     geom_line(lwd=.2,color='firebrick') +
#     scale_y_reverse() +
#     ylab('RT (cycles)') +
#     xlab(paste0('Position on ',slice,' (Mb)')) +
#     coord_cartesian(xlim=c(0,maxX),expand=FALSE) +
#     noLeg + noMarg
#
#   gTop     <- ggarrange(gDSB,gRT,
#                         ncol=1,nrow=2,
#                         align='v')

  return(gX2)
}


#------------------------------------------------------------------------------
#------------------------ FUNCTIONS FOR DISTILL ORIGINS -----------------------
#############################
### getOriginsDF
### KB June 12 2019
### Get the initial origins dataframe
### ARGS:
# origins .tab file (from nextflow)
### OUTPUTS:
# dataframe
# all samplenames
# wt samplenames
getOriginsDFMM10 <- function(fOri){
  oriDF <- fread(fOri,header=TRUE, fill=TRUE)
  sampleNames <- gsub("LW_",replacement = "",names(oriDF)[grep("LW",names(oriDF))])

  # in case we have non-WT samples
  #wtNames <- sampleNames[grep("^wt",sampleNames)]
  ## MANUALLY SPECIFY FOR NOW
  wtNames <- paste0('WT_Rep',1:3)

  # set counters
  oriDF$oneOK <- FALSE
  oriDF$allOK <- TRUE
  oriDF$nOK   <- 0

  oriDF$oneWT <- FALSE
  oriDF$allWT <- TRUE
  oriDF$nWT <- 0

  # Check FR and Ori Presence / Absence for all
  for (s in sampleNames){

    oriDF[[paste0("LW_",s)]] <- oriDF[[paste0("LW_",s)]]+1
    oriDF[[paste0("RW_",s)]] <- oriDF[[paste0("RW_",s)]]+1
    oriDF[[paste0("LC_",s)]] <- oriDF[[paste0("LC_",s)]]+1
    oriDF[[paste0("RC_",s)]] <- oriDF[[paste0("RC_",s)]]+1

    oriDF[[paste0("L_",s)]]  <- oriDF[[paste0("LC_",s)]] / (oriDF[[paste0("LW_",s)]] + oriDF[[paste0("LC_",s)]])
    oriDF[[paste0("R_",s)]]  <- oriDF[[paste0("RW_",s)]] / (oriDF[[paste0("RW_",s)]] + oriDF[[paste0("RC_",s)]])

    oriDF[[paste0("LROK_",s)]]  <- oriDF[[paste0("L_",s)]] < 0.5 & oriDF[[paste0("R_",s)]] < 0.5

    oriDF$oneOK              <- oriDF$oneOK | oriDF[[paste0("LROK_",s)]]
    oriDF$allOK              <- oriDF$allOK & oriDF[[paste0("LROK_",s)]]
    oriDF$nOK                <- oriDF$nOK + oriDF[[paste0("LROK_",s)]]

    oriDF$allWT              <- oriDF$allWT & oriDF[[paste0("ori_",s)]]
    oriDF$oneWT              <- oriDF$oneWT | oriDF[[paste0("ori_",s)]]
    oriDF$nWT                <- oriDF$nWT + (oriDF[[paste0("ori_",s)]]>0)
  }

  ## Define ori clusters as > 7Kb
  oriDF$big                                 <- (oriDF$to-oriDF$from)>7000
  oriDF$oriType                             <- 'Individual'
  oriDF$oriType[(oriDF$to-oriDF$from)>7000] <- 'Cluster'

  oriDF$oriType <- factor(oriDF$oriType,levels=c('Individual','Cluster'))

  return(list(ori=oriDF,
              sampleNames=sampleNames,
              wtSamples=wtNames))
}

#############################
### parseOrigins
### KB June 12 2019
### Parse and compare origins
### ARGS:
# origins data frame (from getOriginsDF)
# type = union or hiconf
### OUTPUTS:
# dataframe of select origins
parseOriginsMM10 <- function(dfOri,oType='union'){

  if (oType == 'union'){
    dfOri$wtOK <- ((dfOri$ori_WT_Rep1 + dfOri$ori_WT_Rep2 + dfOri$ori_WT_Rep3) > 0)
    sampleNames <- c('sperm_Rep1',paste0('WT_Rep',1:3))
  }

  if (oType == 'hiconf'){
    dfOri$wtOK <- ((dfOri$ori_WT_Rep1>0)+(dfOri$ori_WT_Rep2>0)+(dfOri$ori_WT_Rep3>0)) > 1
    sampleNames <- c('sperm_Rep1',paste0('WT_Rep',1:3))
  }

  ## Require OK asymmetry in one sample
  dfOri$LRone <- (dfOri$LROK_WT_Rep1 | dfOri$LROK_WT_Rep2 | dfOri$LROK_WT_Rep3)

  selectDF   <- dfOri[dfOri$wtOK & (dfOri$LRone | dfOri$big),]

  selectDF$minL <- pmin(selectDF$L_WT_Rep1,selectDF$L_WT_Rep2,selectDF$L_WT_Rep3)
  selectDF$minR <- pmin(selectDF$R_WT_Rep1,selectDF$R_WT_Rep2,selectDF$R_WT_Rep3)

  oriName <- data.frame(ori=paste(selectDF$cs,
                                  selectDF$from,
                                  selectDF$to,
                                  sep=":"))


  mylist <- vector(mode="list", length=length(sampleNames))
  names(mylist) <- sampleNames

  for (s in sampleNames){
    mylist[[s]] <- oriName$ori[selectDF[[paste0("ori_",s)]]>0]
  }

  ## Make output files
  originsDF <- selectDF
  originsDF$strength <- (toTPM(originsDF$str_WT_Rep1)+
                           toTPM(originsDF$str_WT_Rep2)+
                           toTPM(originsDF$str_WT_Rep3))/3000000*100

  save(originsDF,file=paste0(oType,'_origins.mm10.Rdata'))

  oriAll2file <- originsDF[,c('cs','from','to','strength')]
  fwrite(x = oriAll2file, file = paste0(oType,'_origins.mm10.bedgraph'),quote = FALSE,sep = "\t",col.names = FALSE)
  fwrite(x = originsDF, file = paste0(oType,'_origins.mm10.details.tab'),quote = FALSE,sep = "\t",col.names = TRUE)

  return(list(ori=originsDF))
}

# Like it says ... parse testis v esc origins
parseOriginsTestis_V_ESCs <- function(dfOri,sample='WT_Rep1'){
  
  sampleNames <- c(paste0('WT_Rep',1:3),'ESC_Rep1','ESC_Rep2')
  
  dfOK <- dfOri %>% 
    mutate(isOri  = ori_WT_Rep1 + ori_WT_Rep2 + ori_WT_Rep3 + ori_ESC_Rep1 + ori_ESC_Rep2 > 0,
           LRok   = LROK_WT_Rep1 | LROK_WT_Rep2 | LROK_WT_Rep3 | LROK_ESC_Rep1| LROK_ESC_Rep2,
           name   = paste(cs,from,to,sep=":"),
           strMei = (toTPM(str_WT_Rep1) + toTPM(str_WT_Rep2) + toTPM(str_WT_Rep3))/3000000*100,
           strESC = toTPM(str_ESC_Rep1)/1000000*100) %>%
    dplyr:::filter(isOri & (LRok | big))
  
  ## Make output files
  save(dfOK,file=paste0('ES_v_Mei_origins.mm10.Rdata'))
  
  oriAll2file <- dfOK[,c('cs','from','to','strMei','strESC')]
  fwrite(x = oriAll2file, file = paste0('ES_v_Mei_origins.mm10.bedgraph'),quote = FALSE,sep = "\t",col.names = FALSE)
  
  return(ori=dfOK)
}

######### PLOT SCATTERPLOT & CC FUNCTION
plotOriCC <- function(df,s1,s2,nm1=NULL,nm2=NULL,chk=TRUE,withDensity=FALSE,pointSz=0.1,pointShape=1,pointAlpha=.1){

  colSet1 <- "#AA4499"
  colSet2 <- "#117733"
  
  if (is.null(nm1)){
    nm1<-s1
  }

  if (is.null(nm2)){
    nm2<-s2
  }

  if (chk){
    df <- df[df[[paste0('ori_',s1)]]>0 | df[[paste0('ori_',s2)]]>0,]
  }

  lblN      <- paste0("\nN = ",format(dim(df)[1], big.mark=",", scientific=FALSE))

  lblCC     <- paste("R^2 == ",round(cor(df[[paste0('str_',s1)]],df[[paste0('str_',s2)]], method='spearman')^2,2),sep='')

  df$type <- 'Neither'
  df$type[df[[paste0('ori_',s1)]]>0 & df[[paste0('ori_',s2)]]<=0] <- paste0(nm1,' only')
  df$type[df[[paste0('ori_',s1)]]<=0 & df[[paste0('ori_',s2)]]>0] <- paste0(nm2,' only')
  df$type[df[[paste0('ori_',s1)]]>0 & df[[paste0('ori_',s2)]]>0]  <- 'Both'

  df$type <- factor(df$type,levels <- c('Both', paste0(nm1,' only'), paste0(nm2,' only'), 'Neither'))

  g <- ggplot(df,aes(color=factor(type),fill=factor(type))) +
    geom_point(aes_string(x=paste0('str_',s1),
                          y=paste0('str_',s2)),
               size=pointSz,
               shape=pointShape,
               alpha=pointAlpha) +
    scale_x_log10(labels=fancy_scientific) +
    scale_y_log10(labels=fancy_scientific) +
    annotate(geom='text',
             x=1,
             y=5*10^3,
             hjust=0,
             vjust=0,
             label=lblN,
             size=7*5/14,
             parse=FALSE) +
    annotate(geom='text',
             x=1,
             y=1*10^3,
             hjust=0,
             vjust=0,
             label=lblCC,
             size=7*5/14,
             parse=TRUE) +
    theme(legend.position='none',
          legend.justification=c(0,1),
          legend.title=element_blank(),
          legend.background=element_blank(),
          legend.margin=unit(c(0,0,0,0),'cm')) +
    annotation_logticks(sides='bl',
                        long=unit(0.1,'cm'),
                        mid=unit(0.05,'cm'),
                        short=unit(0.05,'cm'),
                        lwd=.2,
                        size=0.2) +
    scale_fill_manual(values=c('grey50',colSet1,colSet2)) +
    scale_color_manual(values=c('grey30',colSet1,colSet2)) +
    xlab(nm1) +
    ylab(nm2) + 
    coord_cartesian(xlim=c(.6,2*10^4),ylim=c(.6,2*10^4),expand=FALSE)
    #xlab(paste0('Origin efficiency\n(', nm1,')')) +
    #ylab(paste0('Origin efficiency\n(', nm2,')'))

  if (withDensity){
    g <- g + geom_density_2d(aes_string(x=paste0('str_',s1),
                                        y=paste0('str_',s2)),
                             lwd=.2,
                             alpha=1)
  }

  return(g)
}

######### DRAW SCATTERPLOTS
plotScattersOneSetMM10 <- function(dfOri){
  ## Prep Scatters
  cc12Main <- plotOriCC(dfOri,'WT_Rep1','WT_Rep2', 'Origin efficiency (Testis Rep 1)', 'Origin efficiency (Testis Rep 2)',withDensity = FALSE, pointSz = 0.9, pointShape=21, pointAlpha=.1)
  cc12all  <- plotOriCC(dfOri,'WT_Rep1','WT_Rep2', 'WT Rep 1', 'WT Rep 2')
  cc13all  <- plotOriCC(dfOri,'WT_Rep1','WT_Rep3', 'WT Rep 1', 'WT Rep 3')
  cc32all  <- plotOriCC(dfOri,'WT_Rep3','WT_Rep2', 'WT Rep 3', 'WT Rep 2')

  return(list(gMain = cc12Main,
              g12 = cc12all,
              g13 = cc13all,
              g32 = cc32all))
}

######### DRAW SCATTERPLOTS
plotScattersFigMM10 <- function(dfOriU,dfOriH){
  gUS <- plotScattersOneSetMM10(dfOriU)
  gHS <- plotScattersOneSetMM10(dfOriH)

  if (sum(names(gUS) == 'g42') == 0){
    noYtxt <- theme(axis.text.y=element_blank())
    noXtxt <- theme(axis.text.x=element_blank())
    
    g1 <- ggarrange(gUS$g12 + noXtxt, 
                    gUS$g13 + noYtxt + noXtxt, 
                    gUS$g32 + noYtxt + noXtxt, 
                    ncol=3,nrow=1)
    
    gV1 <- ggarrange(gUS$g12 + noXtxt + theme(plot.margin=unit(c(.5,0,0,0),'cm')), 
                     gUS$g13 + noXtxt, 
                     gUS$g32 , 
                     align='v',
                     heights=c(1.1,1,1),
                     ncol=1,nrow=3)
    
    g2 <- ggarrange(gHS$g12, 
                    gHS$g13 + noYtxt, 
                    gHS$g32 + noYtxt, 
                    ncol=3,nrow=1)
    
    gV2 <- ggarrange(gHS$g12 + noXtxt + noYtxt + theme(plot.margin=unit(c(.5,0,0,0),'cm')), 
                     gHS$g13 + noYtxt + noXtxt, 
                     gHS$g32 + noYtxt, 
                     align='v',
                     heights=c(1.1,1,1),
                     ncol=1,nrow=3)

    gScatters <- ggarrange(grid.arrange(g1,top=textGrob('Origins in ANY sample', gp=gpar(fontsize=7))),
                          grid.arrange(g2,top=textGrob('Origins in at least two samples', gp=gpar(fontsize=7))),
                          nrow=2,ncol=1)
    
    gVerticalScatters <- ggarrange(gV1,gV2,
                                   hjust=0,
                                   labels=c('Origins in ANY','Origins in >= 2 samples'),
                                   nrow=1,ncol=2,
                                   font.label = list(size = 8,
                                                     face='bold'))
    
    ggsave(plot = gScatters,
           getIMGname(fname = 'oriRep_Scatters',saveLocation = imgdir, type = 'PNG'),
           width=(3*2.5),height=(2*2.5),
           dpi = 300, units='in')

    ggsave(plot = gScatters,
           getIMGname(fname = 'oriRep_Scatters',saveLocation = imgdir, type = 'PDF'),
           width=(3*2.5),height=(2*2.5),
           dpi = 300, units='in')
  }
  
  return(list(final=gScatters,
              gVertical=gVerticalScatters,
              gV1=gV1,
              gV2=gV2,
              g1=g1,
              g2=g2,
              gMain=gHS$gMain,
              u12=gUS$g12,
              u13=gUS$g13,
              u32=gUS$g32,
              h12=gHS$g12,
              h13=gHS$g13,
              h32=gHS$g32))
}

######### DRAW SCATTERPLOTS TESTIS V ESCs
plotScatters_Testis_V_ESCs <- function(dfOri){
  ## Prep Scatters
  gCC11 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep1 | dfOri$ori_WT_Rep1],'ESC_Rep1','WT_Rep1', 'ESC Rep 1', 'WT Rep 1')
  gCC12 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep1 | dfOri$ori_WT_Rep2],'ESC_Rep1','WT_Rep2', 'ESC Rep 1', 'WT Rep 2')
  gCC13 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep1 | dfOri$ori_WT_Rep3],'ESC_Rep1','WT_Rep3', 'ESC Rep 1', 'WT Rep 3')
  gEEu  <- plotOriCC(dfOri[dfOri$ori_ESC_Rep1 | dfOri$ori_ESC_Rep2],'ESC_Rep1','ESC_Rep2', 'ESC Rep 1', 'ESC Rep 2')
  gEEi  <- plotOriCC(dfOri[dfOri$ori_ESC_Rep1 & dfOri$ori_ESC_Rep2],'ESC_Rep1','ESC_Rep2', 'ESC Rep 1', 'ESC Rep 2')

  gCC21 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep2 | dfOri$ori_WT_Rep1],'ESC_Rep2','WT_Rep1', 'ESC Rep 2', 'WT Rep 1')
  gCC22 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep2 | dfOri$ori_WT_Rep2],'ESC_Rep2','WT_Rep2', 'ESC Rep 2', 'WT Rep 2')
  gCC23 <- plotOriCC(dfOri[dfOri$ori_ESC_Rep2 | dfOri$ori_WT_Rep3],'ESC_Rep2','WT_Rep3', 'ESC Rep 2', 'WT Rep 3')
  
    
  return(list(g1      = gCC11,
              g2      = gCC12,
              g3      = gCC13,
              gESC2WT = gCC21,
              gEEu    = gEEu,
              gEEi    = gEEi,
              gAll    = ggarrange(gCC11,gEEu,gEEi,ncol=1,nrow=3),
              gVWTc   = ggarrange(gCC11 + theme(axis.text=element_blank(),plot.margin=unit(c(.5,0,0,0),'cm')),
                                gCC12 + theme(axis.text=element_blank()),
                                gCC13 + theme(axis.text.y=element_blank()),
                                ncol=1,nrow=3,
                                align='v',
                                heights=c(1.1,1,1)),
              gAllc   = ggarrange(gCC21 + theme(axis.text=element_blank(),plot.margin=unit(c(.5,0,0,0),'cm')),
                                 gEEu + theme(axis.text=element_blank()),
                                 gEEi + theme(axis.text.y=element_blank()),
                                 ncol=1,nrow=3,
                                 align='v',
                                 heights=c(1.1,1,1))))
}

######### PLOT WATSON V CRICK ASYMMETRY AT ORIS
plotWatsonCrickMM10 <- function(dfOri){
  gMin <- ggplot(dfOri,aes(x=minL,y=minR)) + geom_density_2d(lwd=.2) + geom_hline(color='grey60',lty='dashed',yintercept = 0.5,lwd=.2) + geom_vline(color='grey60',lty='dashed',xintercept = 0.5,lwd=.2) + geom_point(color='grey50',size=.05,alpha=.1)  + facet_grid(.~oriType) + theme(legend.position='none',strip.background=element_blank(),strip.text=element_blank()) + geom_text(aes(label=oriType),x=0.1,y=1,size=7*5/14,check_overlap = TRUE) + coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + xlab('Left of center asymmetry\n(Crick/Watson+Crick)')+ ylab('Right of center asymmetry\n(Watson/Watson+Crick)') + scale_x_continuous(breaks=c(0,.5,1))+ scale_y_continuous(breaks=c(0,.5,1)) + ggtitle('Best value')
  gS1  <- ggplot(dfOri,aes(x=L_WT_Rep1,y=R_WT_Rep1)) + geom_density_2d(lwd=.2) + geom_hline(color='grey60',lty='dashed',yintercept = 0.5,lwd=.2) + geom_vline(color='grey60',lty='dashed',xintercept = 0.5,lwd=.2) + geom_point(color='grey50',size=.05,alpha=.1)  + facet_grid(.~oriType) + theme(legend.position='none',strip.background=element_blank(),strip.text=element_blank()) + geom_text(aes(label=oriType),x=0.1,y=1,size=7*5/14,check_overlap = TRUE) + coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + xlab('Left of center asymmetry\n(Crick/Watson+Crick)')+ ylab('Right of center asymmetry\n(Watson/Watson+Crick)') + scale_x_continuous(breaks=c(0,.5,1))+ scale_y_continuous(breaks=c(0,.5,1))+ ggtitle('Sample 1')
  gS2  <- ggplot(dfOri,aes(x=L_WT_Rep2,y=R_WT_Rep2)) + geom_density_2d(lwd=.2) + geom_hline(color='grey60',lty='dashed',yintercept = 0.5,lwd=.2) + geom_vline(color='grey60',lty='dashed',xintercept = 0.5,lwd=.2) + geom_point(color='grey50',size=.05,alpha=.1)  + facet_grid(.~oriType) + theme(legend.position='none',strip.background=element_blank(),strip.text=element_blank()) + geom_text(aes(label=oriType),x=0.1,y=1,size=7*5/14,check_overlap = TRUE) + coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + xlab('Left of center asymmetry\n(Crick/Watson+Crick)')+ ylab('Right of center asymmetry\n(Watson/Watson+Crick)') + scale_x_continuous(breaks=c(0,.5,1))+ scale_y_continuous(breaks=c(0,.5,1))+ ggtitle('Sample 2')
  gS3  <- ggplot(dfOri,aes(x=L_WT_Rep3,y=R_WT_Rep3)) + geom_density_2d(lwd=.2) + geom_hline(color='grey60',lty='dashed',yintercept = 0.5,lwd=.2) + geom_vline(color='grey60',lty='dashed',xintercept = 0.5,lwd=.2) + geom_point(color='grey50',size=.05,alpha=.1)  + facet_grid(.~oriType) + theme(legend.position='none',strip.background=element_blank(),strip.text=element_blank()) + geom_text(aes(label=oriType),x=0.1,y=1,size=7*5/14,check_overlap = TRUE) + coord_cartesian(xlim=c(0,1),ylim=c(0,1)) + xlab('Left of center asymmetry\n(Crick/Watson+Crick)')+ ylab('Right of center asymmetry\n(Watson/Watson+Crick)') + scale_x_continuous(breaks=c(0,.5,1))+ scale_y_continuous(breaks=c(0,.5,1))+ ggtitle('Sample 3')

  gX4 <- ggarrange(gMin,gS1,gS2,gS3,ncol=1,nrow=4)

  ggsave(plot = gX4,
         getIMGname(fname = 'origin_AsymmetryPlots',saveLocation = imgdir, type = 'PNG'),
         width=4,height=8,
         dpi = 300, units='in')
  ggsave(plot = gX4,
         getIMGname(fname = 'origin_AsymmetryPlots',saveLocation = imgdir, type = 'PDF'),
         width=4,height=8,
         dpi = 300, units='in')

  return(gX4)
}

######### PLOT DINUCLEOTIDES _V ORIGINS
plotDiNT_v_OrisMM10 <- function(dfIn){
  ## ONLY NARROW ORIGINS IN REP1 ON AUTOSOMES
  diNTDF <- as.data.frame(dfIn[dfIn$to-dfIn$from < 7000 &
                                 dfIn$ori_WT_Rep1 &
                                 dfIn$cs %in% paste0('chr',c(1:19)),])

  ## PLOT RATIOS OF DINTs FROM FLANK TO CENTER
  diNTs <- c('GC','CpC','CpG','GpG','GpC','CpA','CpT','GpA','GpT','ApC','ApG','TpC','TpG','ApA','ApT','TpA','TpT')
  for (diNT in diNTs){
    diNTDF[[paste0(diNT,'norm')]]        <- diNTDF[[diNT]]/(diNTDF$to-diNTDF$from)
    diNTDF[[paste0(diNT,'500norm')]]     <- diNTDF[[paste0(diNT,'500')]]/500
    diNTDF[[paste0(diNT,'centerRatio')]] <- log2(diNTDF[[paste0(diNT,'500norm')]]/diNTDF[[paste0(diNT,'norm')]])
  }

  m <- reshape2:::melt.data.frame(diNTDF,id.vars = c('cs','from','to'),
                                  measure.vars = paste0(diNTs,'centerRatio'))

  cgOrder <- c('GC%', 'CpC','CpG','GpG','GpC','CpA','CpT','GpA','GpT','ApC','ApG','TpC','TpG','ApA','ApT','TpA','TpT')

  m$variable <- gsub(m$variable,pattern = 'centerRatio',replacement = '')
  m$variable <- gsub(m$variable,pattern = 'GC',replacement = 'GC%')

  m$variable <- factor(m$variable,levels=cgOrder)

  gDiNT_CenterRatio <- ggplot(m,aes(x=variable,y=value)) +
    geom_hline(yintercept=0,lwd=.2) +
    geom_violin(lwd=.2,fill='grey80',adjust=2,scale = 'width',fill='grey99',trim = TRUE) +
    geom_boxplot(width=.15,lwd=.2,color='black',fill='pink',outlier.size=.1,outlier.alpha=.1,notch=TRUE) +
    theme(legend.position='none') +
    ylab('Central enrichment\nlog2 (center / flanks)') + xlab('') +
    theme(legend.position='none',axis.text.x = element_text(angle = 90)) +
    geom_vline(xintercept=1.5,lty='dashed',lwd=.2)


  ## PLOT CCs V STRENGTH
  cType <- "500norm"
  nt <- c('G','A','T','C')

  ntCol <- c(paste0('Gp',nt,cType),paste0('Ap',nt,cType),paste0('Tp',nt,cType),paste0('Cp',nt,cType))

  mCC <- cor(diNTDF[,c('str_WT_Rep1','str_WT_Rep2','str_WT_Rep3','AT500','GC500',ntCol)],method='spearman')

  d <- as.data.frame(mCC)[,1:3]
  d$type <- rownames(d)

  mD <- reshape2:::melt.data.frame(d,measure.vars = names(d)[1:3],id.vars = c('type'))

  mD <- mD[mD$variable != 'str_WT_Rep2',]
  mD <- mD[mD$variable != 'str_WT_Rep3',]
  mD <- mD[mD$type     != 'str_WT_Rep1',]
  mD <- mD[mD$type     != 'str_WT_Rep2',]
  mD <- mD[mD$type     != 'str_WT_Rep3',]
  mD <- mD[mD$type     != 'AT500',]

  mD$type <- gsub(mD$type,pattern = cType,replacement = '')
  mD$type <- gsub(mD$type,pattern = '500',replacement = '%')

  cType <- '';
  cgOrder <- c('GC%','CpC','CpG','GpG','GpC','CpA','CpT','GpA','GpT','ApC','ApG','TpC','TpG','ApA','ApT','TpA','TpT')
  mD$type <- factor(mD$type,levels=cgOrder)

  gDiNTCC <- ggplot(mD,aes(x=type,y=value,fill=value)) +
    geom_bar(color='black',lwd=.2,alpha=.7,stat='identity',width=.75) +
    geom_vline(xintercept=1.5,lty='dashed',lwd=.2) +
    geom_hline(yintercept=0,lwd=.2) +
    coord_cartesian(ylim=c(-1,1)) +
    ylab('Correlation with origin strength\n(Spearman R)') +
    xlab('') +
    scale_fill_gradient2(low='forestgreen',midpoint=0,mid='white',high='firebrick') +
    theme(legend.position='none',axis.text.x = element_text(angle = 90,vjust=.5,hjust=1))

  gDiNTplot <- ggarrange(gDiNT_CenterRatio + theme(axis.text.x = element_blank()),
                         gDiNTCC,
                         nrow=2,ncol=1,
                         align='v',
                         heights=c(4,5))

  #save(selectDF,file='oriRepNarrow.Rdata')

  ######################## PLOT THE COMPOSITE FIGURE
  ggsave(plot = gDiNTplot,
         getIMGname(fname = 'origins_V_diNTs',saveLocation = imgdir, type = 'PNG'),
         width=3.5,height=3,
         dpi = 300, units='in')
  ggsave(plot = gDiNTplot,
         getIMGname(fname = 'origins_V_diNTs',saveLocation = imgdir, type = 'PDF'),
         width=3.5,height=3,
         dpi = 300, units='in')

  return(list(both=gDiNTplot,
              gCenter=gDiNT_CenterRatio,
              gCC=gDiNTCC))
}

######### PLOT RAW ORI COUNTS PER SAMPLE
plotOriCounts <- function(dfOri,cSamples=NULL){

  if (is.null(cSamples)){
    cSamples <- grep(x=names(dfOri),pattern = 'ori_')
  }

  sCols <- paste0('ori_',cSamples)

  print(c('cs','from','to',sCols))

  mCount  <- reshape2:::melt.data.frame(data = as.data.frame(dfOri)[,c('cs','from','to',sCols)],
                                        measure.vars=sCols,
                                        id.vars=c('cs','from','to'))

  mCount$sample <- gsub(x= mCount$variable,pattern = 'ori_',replacement = '')

  countDF <- as.data.frame(aggregate(mCount$value,
                                     by=list(sample=mCount$sample),
                                     FUN=function(x){sum(x>0)}))

  gOriCounts <- ggplot(countDF,aes(x=sample,y=x,
                                   label=paste0('N = ',format(x = x, big.mark = ",", scientific=FALSE)))) +
    geom_bar(stat='identity') +
    geom_text(y=100,color='white',hjust=0,size=8*5/14) +
    coord_flip() + xlab('') + ylab('Number of detected origins')

  ggsave(plot = gOriCounts,
         getIMGname(fname = 'oriCounts',saveLocation = imgdir, type = 'PNG'),
         width=3,height=2,
         dpi = 300, units='in')
  ggsave(plot = gOriCounts,
         getIMGname(fname = 'oriCounts',saveLocation = imgdir, type = 'PDF'),
         width=3,height=2,
         dpi = 300, units='in')

  return(gOriCounts)
}

######### PARSE ORIGINS FOR COMPARISON WITH RNA hydrolysis CTRL
parseWT4Oris <- function(dfOri){

  dfOri$wtOK <- ((dfOri$ori_WT_Rep4>0)+(dfOri$ori_WT_Rep2>0)+(dfOri$ori_WT_Rep3>0)) > 1

  dfOri$LRwt <- (dfOri$LROK_WT_Rep1 | dfOri$LROK_WT_Rep2 | dfOri$LROK_WT_Rep3)

  selectDF   <- dfOri[(dfOri$ori_WT_Rep4 > 0) &
                        (dfOri$ori_WT_Rep1 > 0 | dfOri$ori_WT_Rep2 > 0 | dfOri$ori_WT_Rep3 > 0) &
                        (dfOri$LROK_WT_Rep4 | dfOri$LRwt | dfOri$big),]

  oriName <- data.frame(ori=paste(selectDF$cs,
                                  selectDF$from,
                                  selectDF$to,
                                  sep=":"))

  oType <- 'WT4'

  sampleNames <- paste0('WT_Rep',1:4)

  mylist <- vector(mode="list", length=length(sampleNames))
  names(mylist) <- sampleNames

  for (s in sampleNames){
    mylist[[s]] <- oriName$ori[selectDF[[paste0("ori_",s)]]>0]
  }

  ## Show Ori Overlaps
  graphics.off()
  png(getIMGname(fname = paste0(oType,'_origins_UpSet'),saveLocation = imgdir, type = 'PNG'),
      width=7,height=6, res = 300, units='in')
  upset(fromList(mylist), order.by = "freq")
  dev.off()

  vennSVG <- getIMGname(fname = paste0(oType,'_origins_Venn'),saveLocation = imgdir, type = 'SVG')
  venn <- plotVenn(mylist[paste0('WT_Rep',1:3)],
                   nCycles=5000,
                   systemShow=FALSE,
                   fontScale=2.25,
                   setColors=c('orange','grey','firebrick'),
                   borderWidth=3)
  showSVG(nVennObj = venn, outFile = vennSVG)

  ## Make output files
  originsDF <- selectDF
  originsDF$strength <- (toTPM(originsDF$str_WT_Rep1)+
                           toTPM(originsDF$str_WT_Rep2)+
                           toTPM(originsDF$str_WT_Rep3))/3000000*100

  save(originsDF,file=paste0(oType,'_origins.mm10.Rdata'))

  oriAll2file <- originsDF[,c('cs','from','to','strength')]
  fwrite(x = oriAll2file, file = paste0(oType,'_origins.mm10.bedgraph'),quote = FALSE,sep = "\t",col.names = FALSE)

  return(originsDF)
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
plotOriOverlaps <- function(dfOri,sampleNames=NULL){
  library(grImport2)
  library(ggplotify)

  oriName <- data.frame(ori=paste(dfOri$cs,
                                  dfOri$from,
                                  dfOri$to,
                                  sep=":"))

  mylist <- vector(mode="list", length=length(sampleNames))
  names(mylist) <- gsub(x=sampleNames,pattern = 'wt_', replacement = '')

  for (s in sampleNames){
    listNM           <- gsub(x=s, pattern = 'wt_', replacement = '')
    mylist[[listNM]] <- oriName$ori[dfOri[[paste0("ori_",s)]]>0]
  }

  ## Show Ori Overlaps
  ggg <- as.grob(upset(fromList(mylist),
                       order.by = "freq",
                       matrix.color = 'red',
                       mainbar.y.label = 'Origins (count)',
                       sets.x.label = 'Origins per sample',
                       shade.color = 'pink'))

  gCounts    <- plotOriCounts(dfOri,sampleNames)
  gUpSetMain <- ggg$children[[1]]$grobs[[1]]
  gUpSetInt  <- ggg$children[[1]]$grobs[[2]]
  gUpSetTot  <- ggg$children[[2]]$grobs[[1]]

  noMarg     <- theme(plot.margin = unit(c(0,0,0,0),'cm'))

  ggUpSet    <- ggarrange(ggplot()  ,
                          gUpSetMain ,
                          gUpSetTot      ,
                          gUpSetInt  ,
                          ncol=2,nrow=2,
                          widths=c(1,3),
                          heights=c(3,1),
                          align='hv')

  return(ggUpSet)
}

######################################################################
## Compare features relative to RT
plotFeaturesVRT <- function(x,field2use,RT,nQ,rev=TRUE,eLbl=NULL,nR=10){

  if (is.null(eLbl)){
    eLbl <- field2use
  }

  dfMn <- aggregate(x[[field2use]],
                    by=list(RT=convertToQuantiles(x[[RT]],
                                                  nQ = nQ,
                                                  revLbl=rev,
                                                  numericQuantiles = TRUE)),
                    FUN=mean)

  dfLen <- aggregate(x[[field2use]],
                     by=list(RT=convertToQuantiles(x[[RT]],
                                                   nQ = nQ,
                                                   revLbl=rev,
                                                   numericQuantiles = TRUE)),
                     FUN=function(x){sum(x>0)})
  names(dfLen) <- c('RT','N')

  df <- plyr:::join(dfMn,dfLen,by='RT')

  df$xRaw  <- df$x
  df$x     <- standardizeMNSD(df$x)

  df$type  <- 'RT'

  if (nR > 0){
    for (n in 1:nR){
      dfRiMn <- aggregate(sample(x[[field2use]]),
                          by=list(RT=convertToQuantiles(x[[RT]],
                                                        nQ = nQ,
                                                        revLbl=TRUE,
                                                        numericQuantiles = TRUE)),
                          FUN=mean)

      dfRiLen <- aggregate(sample(x[[field2use]]),
                           by=list(RT=convertToQuantiles(x[[RT]],
                                                         nQ = nQ,
                                                         revLbl=TRUE,
                                                         numericQuantiles = TRUE)),
                           FUN=function(x){sum(x>0)})
      names(dfRiLen) <- c('RT','N')

      dfRi <- plyr:::join(dfRiMn,dfRiLen,by='RT')
      dfRi$type  <- 'random'

      dfRi$xRaw  <- dfRi$x
      dfRi$x     <- standardizeMNSD(dfRi$x)

      if (n == 1){
        dfR <- dfRi
      }else{
        dfR <- rbind(dfR,dfRi)
      }
    }
    df <- rbind(df,dfR)
  }



  df$type <- factor(df$type,levels=c('RT','random'))
  df$expt <- eLbl

  gRet <- ggplot(df,aes(x=RT,y=x,color=type,fill=type)) +
    scale_color_manual(values=c('orange','grey80')) +
    scale_fill_manual(values=c('orange','grey90')) +
    scale_size_manual(values=c(0.6,0.1)) +
    geom_point(aes(size=type)) +
    geom_smooth(lwd=.2) +
    theme(legend.position='none') +
    ylab(eLbl)
  #print(df)
  return(list(fig=gRet,
              data=df))
}

######################################################################
## Compare features relative to RT (multiple feats)
plotFeaturesVRTMulti <- function(iDF,
                                 iFields,
                                 names=NULL,
                                 iRT,
                                 q,
                                 nRandom=10,
                                 seAlpha=.4,
                                 smoothSpan=2,
                                 pointSz=.3,
                                 pCol='orange',
                                 ylbl=NULL,
                                 rev=FALSE,
                                 showDiagonal=FALSE,
                                 showZero=TRUE,
                                 colorByType=TRUE){

  if (exists('allDataFPlot')){rm('allDataFPlot')}
  if (exists('iD')){rm('iD')}

  for (f in iFields){
    iD <- plotFeaturesVRT(iDF,f,iRT,q,rev=rev,nR=nRandom)

    if (exists('allDataFPlot')){
      #print(dim(allDataFPlot))
      allDataFPlot <- rbind(allDataFPlot,iD$data)
    }else{
      allDataFPlot <- iD$data
    }
  }

  yMax <- quantile(allDataFPlot$x,.975,na.rm=TRUE)
  yMin <- -1*yMax

  ## Replace field names with labels
  if (is.null(names)){
    allDataFPlot$expt <- factor(allDataFPlot$expt,levels=iFields)
  }else{
    nmDF <- data.frame(expt=iFields,name=names)
    ad   <- plyr:::join(allDataFPlot,nmDF,by='expt')
    ad$expt <- factor(ad$name,levels=names)
    allDataFPlot <- ad
    print(head(ad))
  }

  if (is.null(ylbl)){
    yLabelText <- 'Relative density'
  }else{
    yLabelText <- bquote(.(ylbl)*"\nRelative density")
  }

  if (colorByType){
    types  <- c('affySeq','PRDM9ChIPSeq','h3k4m3','spo11','SSDS','RPA','RPA-SSDS','CO','NCO','IH','other','random')

    colGG <- scale_color_manual(values=c('affySeq'      = 'royalblue1',
                                         'PRDM9ChIPSeq' = 'dodgerblue4',
                                         'h3k4m3'       = 'orange',
                                         'spo11'        = 'magenta',
                                         'SSDS'         = 'forestgreen',
                                         'RPA-SSDS'     = 'brown',
                                         'CO'           = 'palevioletred2',
                                         'NCO'          = 'black',
                                         'IH'           = 'grey50',
                                         'other'        = 'cyan',
                                         'random'       = 'grey80'))

    fillGG <- scale_fill_manual(values=c('affySeq' = 'royalblue1',
                                         'PRDM9ChIPSeq'  = 'dodgerblue4',
                                         'h3k4m3'       = 'orange',
                                         'spo11'        = 'plum1',
                                         'SSDS'         = 'forestgreen',
                                         'RPA-SSDS'     = 'brown',
                                         'CO'           = 'palevioletred2',
                                         'NCO'          = 'black',
                                         'IH'          = 'grey50',
                                         'other'        = 'cyan',
                                         'random'       = 'grey80'))

    sizeGG <- scale_size_manual(values=c('affySeq' = pointSz,
                                         'PRDM9ChIPSeq'  = pointSz,
                                         'h3k4m3'       = pointSz,
                                         'spo11'        = pointSz,
                                         'SSDS'         = pointSz,
                                         'RPA-SSDS'     = pointSz,
                                         'CO'           = pointSz,
                                         'NCO'          = pointSz,
                                         'IH'           = pointSz,
                                         'other'        = pointSz,
                                         'random'       = 0))

    alphaGG <- scale_alpha_manual(values=c('affySeq' = .3,
                                           'PRDM9ChIPSeq'= .3,
                                           'h3k4m3'= .3,
                                           'spo11'= .3,
                                           'SSDShop2' = .3,
                                           'SSDS'= .3,
                                           'RPA-SSDS'= .3,
                                           'CO'= .3,
                                           'NCO'= .3,
                                           'IH'= .3,
                                           'other'= .3,
                                           'random'= 0))

    names(allDataFPlot)[which(names(allDataFPlot) == 'type')]        <- 'itype'

    allDataFPlot$type                                              <- 'other'
    allDataFPlot$type[grep(pattern = 'affySeq',allDataFPlot$expt)]      <- 'affySeq';
    allDataFPlot$type[grep(pattern = 'PRDM9ChIPSeq',allDataFPlot$expt)] <- 'PRDM9ChIPSeq';
    allDataFPlot$type[grep(pattern = 'h3k4m3',allDataFPlot$expt)]       <- 'h3k4m3';
    allDataFPlot$type[grep(pattern = 'spo11',allDataFPlot$expt)]        <- 'spo11';
    allDataFPlot$type[grep(pattern = 'SSDS',allDataFPlot$expt)]         <- 'SSDS';
    allDataFPlot$type[grep(pattern = 'RPA',allDataFPlot$expt)]         <- 'RPA-SSDS';
    allDataFPlot$type[grep(pattern = 'RR',allDataFPlot$expt)]        <- 'CO';
    allDataFPlot$type[grep(pattern = 'Crossover',allDataFPlot$expt)] <- 'CO';
    allDataFPlot$type[grep(pattern = 'CO',allDataFPlot$expt)]        <- 'CO';
    allDataFPlot$type[grep(pattern = 'NCO',allDataFPlot$expt)]       <- 'NCO';
    allDataFPlot$type[grep(pattern = 'CONCOhh',allDataFPlot$expt)]   <- 'IH';
    allDataFPlot$type[grep(pattern = 'IH',allDataFPlot$expt)]        <- 'IH';
    allDataFPlot$type[grep(pattern = 'Non-C',allDataFPlot$expt)]     <- 'NCO';
    allDataFPlot$type[grep(pattern = 'random',allDataFPlot$itype)]   <- 'random';

    allDataFPlot$type <- factor(allDataFPlot$type,levels=types)

    titlez  <- c('affySeq','PRDM9ChIPSeq','H3K4m3','SPO11','DMC1-SSDS',
                 'DMC1-SSDS (wt Rep1)','DMC1-SSDS (wt Rep2)',
                 'DMC1-SSDS (Hop2KO)','DMC1-SSDS (Prdm9KO)',
                 'SSDS','RPA','RPA-SSDS','CO','NCO','IH','other','random')

    allDataFPlot$title                                                   <- 'other'
    allDataFPlot$title[grep(pattern = 'affySeq',allDataFPlot$expt)]      <- 'affySeq';
    allDataFPlot$title[grep(pattern = 'PRDM9ChIPSeq',allDataFPlot$expt)] <- 'PRDM9ChIPSeq';
    allDataFPlot$title[grep(pattern = 'h3k4m3',allDataFPlot$expt)]       <- 'H3K4M3';
    allDataFPlot$title[grep(pattern = 'spo11',allDataFPlot$expt)]        <- 'SPO11-oligos';
    allDataFPlot$title[grep(pattern = 'SSDST1',allDataFPlot$expt)]       <- 'DMC1-SSDS (wt Rep1)';
    allDataFPlot$title[grep(pattern = 'SSDST2',allDataFPlot$expt)]       <- 'DMC1-SSDS (wt Rep2)';
    allDataFPlot$title[grep(pattern = 'SSDShop2',allDataFPlot$expt)]     <- 'DMC1-SSDS (Hop2KO)';
    allDataFPlot$title[grep(pattern = 'SSDSPrKO',allDataFPlot$expt)]     <- 'DMC1-SSDS (Prdm9KO)';
    allDataFPlot$title[grep(pattern = '^SSDS$',allDataFPlot$expt)]       <- 'DMC1-SSDS';
    allDataFPlot$title[grep(pattern = 'RPA',allDataFPlot$expt)]          <- 'RPA-SSDS';
    allDataFPlot$title[grep(pattern = 'RR',allDataFPlot$expt)]        <- 'CO';
    allDataFPlot$title[grep(pattern = 'Crossover',allDataFPlot$expt)] <- 'CO';
    allDataFPlot$title[grep(pattern = 'CO',allDataFPlot$expt)]        <- 'CO';
    allDataFPlot$title[grep(pattern = 'NCO',allDataFPlot$expt)]       <- 'NCO';
    allDataFPlot$title[grep(pattern = 'Non-C',allDataFPlot$expt)]     <- 'NCO';
    allDataFPlot$title[grep(pattern = 'CONCOhh',allDataFPlot$expt)]     <- 'Inter-homolog repair';
    allDataFPlot$title[grep(pattern = 'IH',allDataFPlot$expt)]     <- 'Inter-homolog repair';
    allDataFPlot$title[grep(pattern = 'random',allDataFPlot$ititle)]   <- 'random';

    allDataFPlot$titlez <- factor(allDataFPlot$title,levels=titlez)

  }else{
    allDataFPlot$type  <- factor(allDataFPlot$type,levels = c('random','RT'))
    allDataFPlot$title <- allDataFPlot$type
    colz <- c('grey80',pCol)
    sizez  <- c(0,0.6)
    alphaz <- c(0,1)
    colGG  <- scale_color_manual(c('random' = 'grey80','RT' = pCol))
    fillGG <- scale_fill_manual(c('random' = 'grey80','RT' = pCol))
    sizeGG <- scale_size_manual(c('random' = 0,'RT' = 0.8))
    alphaGG <- scale_alpha_manual(c('random' = 0,'RT' = 1))
  }

  gFig <- ggplot(allDataFPlot,aes(x=RT,y=x,color=type,fill=type)) +
    colGG +
    fillGG +
    sizeGG +
    alphaGG +
    geom_point(aes(size=type,alpha=type)) +
    geom_smooth(span=smoothSpan,alpha=seAlpha,lwd=.3,method='loess') +
    theme(legend.position='none') +
    ylab(yLabelText) +
    facet_wrap(~expt,nrow=1) +
    theme(legend.position='none',
          strip.background=element_blank(),
          strip.text = element_blank()) +
    geom_text(x=q/2,y=Inf,
              hjust=0.5,vjust=1.2,
              size=7*5/14,
              color='black',
              check_overlap = TRUE,
              aes(label=title)) +
    coord_cartesian(ylim=c(yMin,yMax)) +
    scale_x_continuous(breaks=seq(1,q,length.out = 5),
                       labels = c('0','','0.5','','1')) +
    xlab('Relative RT (percentile)')

  gFigNoFit <- ggplot(allDataFPlot,aes(x=RT,y=x,color=type,fill=type)) +
    colGG +
    fillGG +
    sizeGG +
    alphaGG +
    geom_point(aes(size=type,alpha=type)) +
    geom_smooth(data=allDataFPlot[allDataFPlot$type=='random',],span=2,lwd=.3,method='loess') +
    theme(legend.position='none') +
    ylab(yLabelText) +
    facet_wrap(~expt,nrow=1) +
    theme(legend.position='none',
          strip.background=element_blank(),
          strip.text = element_blank()) +
    geom_text(x=q/2,y=Inf,
              hjust=0.5,vjust=1.2,
              size=7*5/14,
              color='black',
              check_overlap = TRUE,
              aes(label=title)) +
    coord_cartesian(ylim=c(yMin,yMax)) +
    scale_x_continuous(breaks=seq(1,q,length.out = 5),
                       labels = c('0','','0.5','','1')) +
    xlab('Relative RT (percentile)')

  gFigRaw <- ggplot(allDataFPlot,aes(x=RT,y=xRaw,color=type,fill=type)) +
    colGG +
    fillGG +
    sizeGG +
    alphaGG +
    geom_point(aes(size=type,alpha=type),shape=21,fill='grey50') +
    geom_smooth(span=2,lwd=.3,method='loess') +
    theme(legend.position='none') +
    ylab('Relative density') +
    facet_wrap(~expt,nrow=1,scales='free_y') +
    theme(legend.position='none',
          strip.background=element_blank(),
          strip.text = element_blank()) +
    geom_text(x=q/2,y=Inf,
              hjust=0.5,vjust=1.2,
              size=7*5/14,
              color='black',
              check_overlap = TRUE,
              aes(label=title)) +
    scale_x_continuous(breaks=seq(1,q,length.out = 5),
                       labels = c('0','','0.5','','1')) +
    xlab('Relative RT (percentile)')

  ## Start of Alt fig ! Keep ;) !
  gFigAlt <- ggplot(allDataFPlot,aes(x=RT,y=x,color=type,fill=type)) +
    colGG +
    fillGG +
    sizeGG +
    alphaGG

  if (showDiagonal){
    #diagXY = 1.5
    diagXY = 2
    if (rev){
      gFigAlt <- gFigAlt + geom_line(data=data.frame(x=c(0,q),y=c(diagXY,-diagXY),type='SSDShh'),
                                    aes(x=x,y=y),
                                    lwd=.2,lty='dashed',color='red')
    }else{
      gFigAlt <- gFigAlt + geom_line(data=data.frame(x=c(0,q),y=c(-diagXY,diagXY),type='SSDShh'),
                                    aes(x=x,y=y),
                                    lwd=.2,lty='dashed',color='red')
    }
  }

  if (showZero){
      gFigAlt <- gFigAlt + geom_line(data=data.frame(x=c(0,q),y=c(0,0),type='SSDShh'),
                                     aes(x=x,y=y),
                                     lwd=.2,lty='dashed',color='grey60')
  }
  
  gFigAlt <- gFigAlt + geom_point(aes(size=type,alpha=type)) +
    geom_smooth(span=smoothSpan,alpha=seAlpha,lwd=.5,method='loess') +
    theme(legend.position='none') +
    ylab(yLabelText) +
    theme(legend.position='none',
          strip.background=element_blank(),
          strip.text = element_blank()) +
    coord_cartesian(ylim=c(yMin,yMax)) +
    scale_x_continuous(breaks=seq(1,q,length.out = 5),
                       labels = c('0','','0.5','','1')) +
    xlab('Replication Timing')

  return(retData=list(fig=gFig,
                      figRaw=gFigRaw,
                      figAlt=gFigAlt,
                      data=allDataFPlot))
}

#### processAndcleanStatsFile
getModelStats <- function(statFile=NULL,bestModelCutoff= 0.9985){

  ## import metrics file
  statDF       <- fread(statFile,header=TRUE)

  ## set species
  statDF$species                                                  <- 'human'
  statDF$species[grep(pattern = 'mm10', as.character(statDF$bg))] <- 'mouse'

  ## Replace bits of ori name
  statDF$ori[statDF$ori == 'ESC_Rep1.peaks.FRfiltered.bedgraph']     <- 'OriSSDS(ESC)'
  statDF$ori[statDF$ori == 'WT_Rep1.peaks.FRfiltered.bedgraph']      <- 'OriSSDS(WT1)'

  statDF$ori <- gsub(pattern = '.mm10.bedgraph',            replacement = '',            as.character(statDF$ori))
  statDF$ori <- gsub(pattern = '.hg38.bedgraph',            replacement = '',            as.character(statDF$ori))
  statDF$ori <- gsub(pattern = '.bedgraph',                 replacement = '',            as.character(statDF$ori))
  statDF$ori <- gsub(pattern = 'randomized_hiconf_origins', replacement = 'OriSSDS(r)',  as.character(statDF$ori))
  statDF$ori <- gsub(pattern = 'hiconf_origins',            replacement = 'OriSSDS(hi)', as.character(statDF$ori))
  statDF$ori <- gsub(pattern = 'union_origins',             replacement = 'OriSSDS(u)',  as.character(statDF$ori))

  ## Rename samples

  statDF$bg[statDF$bg == "2C_ALL_hsP4_2C_NA_REC8.hg38.RT.forModel.bedgraph"]               <- "2C"
  statDF$bg[statDF$bg == "MeiS_ALL_hsS1Ad_2to4C_SCP3_NA.hg38.RT.forModel.bedgraph"]        <- "MeiS(S)a"
  statDF$bg[statDF$bg == "MeiS_ALL_hsS1_2to4C_SCP3_NA.hg38.RT.forModel.bedgraph"]          <- "MeiS(S)b"
  statDF$bg[statDF$bg == "MeiS_ALL_hsS1_2to4C_HORMAD1DMRT1_NA.hg38.RT.forModel.bedgraph"]  <- "MeiS(H/D)"
  statDF$bg[statDF$bg == "MeiS_weakH1_hsS1_2to4C_HORMAD1_DMRT1.hg38.RT.forModel.bedgraph"] <- "MeiS(H1w)"
  statDF$bg[statDF$bg == "RPE1_MID1_RT_hg38.forModel.bedgraph"]                            <- "RPE1a"
  statDF$bg[statDF$bg == "RPE1_MID2_RT_hg38.forModel.bedgraph"]                            <- "RPE1b"
  statDF$bg[statDF$bg == "RPE1_MID3_RT_hg38.forModel.bedgraph"]                            <- "RPE1c"
  statDF$bg[statDF$bg == "RT_BG02_ESC_Int51992745_hg38.forModel.bedgraph"]                 <- "BG02ESCa"
  statDF$bg[statDF$bg == "RT_BG02_ESC_Int87960943_hg38.forModel.bedgraph"]                 <- "BG02ESCb"
  statDF$bg[statDF$bg == "RT_CyT49_Liver_D16_Int81158282_hg38.forModel.bedgraph"]          <- "Liver"
  statDF$bg[statDF$bg == "RT_FM01-154-001_Myoblast_Int58331187_hg38.forModel.bedgraph"]    <- "Myoblast"
  statDF$bg[statDF$bg == "RT_GM06990_Lymphoblastoid_Ext54054609_hg38.forModel.bedgraph"]   <- "Lymphoblastoid"
  statDF$bg[statDF$bg == "RT_H7_ESC_Ext35479608_hg38.forModel.bedgraph"]                   <- "H7 ESC"
  statDF$bg[statDF$bg == "RT_H9_Neural_Progenitor_Int89790558_hg38.forModel.bedgraph"]     <- "H9 NPC"
  statDF$bg[statDF$bg == "RT_IMR90_Lung_Fibroblast_Int49605910_hg38.forModel.bedgraph"]    <- "IMR90 Lung Fibro"
  statDF$bg[statDF$bg == "RT_U2OS_Bone_Epithelial_Int66343918_hg38.forModel.bedgraph"]     <- "U2OS Bone Epi"

  statDF$bg[statDF$bg == "2C_ALL_S2_2C_NA_NA.mm10.RT.forModel.bedgraph"]                   <-  "2C"
  statDF$bg[statDF$bg == "BCell_Tubbs_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]              <-  "Bcell"
  statDF$bg[statDF$bg == "CD8_Yehuda_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]               <-  "CD8"
  statDF$bg[statDF$bg == "E14_Dey_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]                  <-  "ESC(b)"
  statDF$bg[statDF$bg == "ESC_ALL_S1_2to4C_NA_NA.mm10.RT.forModel.bedgraph"]               <-  "ESC(a)"
  statDF$bg[statDF$bg == "Lep_ALL_S1_4C_SCP3STRA8_NA.mm10.RT.forModel.bedgraph"]           <-  "4C"
  statDF$bg[statDF$bg == "MeiS_ALL_S3_2to4C_STRA8_DMRT1.mm10.RT.forModel.bedgraph"]        <-  "MeiS(a)"
  statDF$bg[statDF$bg == "MeiS_ALL_S4_2to4C_STRA8_DMRT1.mm10.RT.forModel.bedgraph"]        <-  "MeiS(b)"
  statDF$bg[statDF$bg == "MeiS_ALL_S5_2to4C_STRA8_DMRT1.mm10.RT.forModel.bedgraph"]        <-  "MeiS(c)"
  statDF$bg[statDF$bg == "MeiS_SPO11_S1_2to4C_STRA8_DMRT1.mm10.RT.forModel.bedgraph"]      <-  "MeiS(SPO11)"
  statDF$bg[statDF$bg == "MeiS_wtLM_S1_2to4C_STRA8_DMRT1.mm10.RT.forModel.bedgraph"]       <-  "MeiS(wtLM)"
  statDF$bg[statDF$bg == "MeiS_BxCEARLY_R1_2to4C_SCP3_NA.mm10.RT.forModel.bedgraph"]       <-  "MeiS(BxC)a"
  statDF$bg[statDF$bg == "MeiS_BxCLATE_R1_2to4C_SCP3_NA.mm10.RT.forModel.bedgraph"]        <-  "MeiS(BxC)b"
  statDF$bg[statDF$bg == "MeiS_BxCMID_R1_2to4C_SCP3_NA.mm10.RT.forModel.bedgraph"]         <-  "MeiS(BxC)c"
  statDF$bg[statDF$bg == "MeiS_EARLY_S2_2to4C_SCP3yH2AX_NA.mm10.RT.forModel.bedgraph"]     <-  "MeiS(T2)a"
  statDF$bg[statDF$bg == "MeiS_EARLY_S3_2to4C_SCP3STRA8_NA.mm10.RT.forModel.bedgraph"]     <-  "MeiS(T2)b"
  statDF$bg[statDF$bg == "MeiS_LATE_S2_2to4C_SCP3yH2AX_NA.mm10.RT.forModel.bedgraph"]      <-  "MeiS(T4)"
  statDF$bg[statDF$bg == "MeiS_MID_S2_2to4C_SCP3yH2AX_NA.mm10.RT.forModel.bedgraph"]       <-  "MeiS(T3)"
  statDF$bg[statDF$bg == "MeiS_VERYEARLY_S1_2to4C_SCP3_NA.mm10.RT.forModel.bedgraph"]      <-  "MeiS(T1)"
  statDF$bg[statDF$bg == "PGC_Yehuda_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]               <-  "PGC"
  statDF$bg[statDF$bg == "PREB_Yehuda_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]              <-  "PRE-B"
  statDF$bg[statDF$bg == "SPGB_ALL_S2_2to4C_STRA8DMRT1hi_NA.mm10.RT.forModel.bedgraph"]    <-  "SpgB(a)"
  statDF$bg[statDF$bg == "SPGB_ALL_S4_2to4C_STRA8DMRT1_NA.mm10.RT.forModel.bedgraph"]      <-  "SpgB(b)"
  statDF$bg[statDF$bg == "SPGI_ALL_S1_2to4C_DMRT1hi_STRA8.mm10.RT.forModel.bedgraph"]      <-  "SpgI(a)"
  statDF$bg[statDF$bg == "SPGI_ALL_S4_2to4C_DMRT1hi_STRA8.mm10.RT.forModel.bedgraph"]      <-  "SpgI(b)"
  statDF$bg[statDF$bg == "SPGI_SPO11_S1_2to4C_DMRT1hi_STRA8.mm10.RT.forModel.bedgraph"]    <-  "SpgI(SPO11)"
  statDF$bg[statDF$bg == "SPGI_wtLM_S1_2to4C_DMRT1hi_STRA8.mm10.RT.forModel.bedgraph"]     <-  "SpgI(wtLM)"
  statDF$bg[statDF$bg == "SPGU_ALL_S2_2to4C_DMRT1_STRA8.mm10.RT.forModel.bedgraph"]        <-  "SpgU(a)"
  statDF$bg[statDF$bg == "SPGU_ALL_S4_2to4C_DMRT1_STRA8.mm10.RT.forModel.bedgraph"]        <-  "SpgU(b)"
  statDF$bg[statDF$bg == "SPGU_SPO11_S1_2to4C_DMRT1_STRA8.mm10.RT.forModel.bedgraph"]      <-  "SpgU(SPO11)"
  statDF$bg[statDF$bg == "SPGU_wtLM_S1_2to4C_DMRT1_STRA8.mm10.RT.forModel.bedgraph"]       <-  "SpgU(wtLM)"
  statDF$bg[statDF$bg == "SSC_Yehuda_R1_NA_NA_NA.mm10.RT.forModel.bedgraph"]               <-  "SSC"
  statDF$bg[statDF$bg == "Sertoli_ALL_S2_2C_GATA_NA.mm10.RT.forModel.bedgraph"]            <-  "Sertoli"
  statDF$bg[statDF$bg == "RT_C127_Mammary.mm10.forModel.bedgraph"]                         <-  "Mammary"
  statDF$bg[statDF$bg == "RT_D3_NPC_EBM9_Int88652090_mm10.forModel.bedgraph"]              <-  "NPC"
  statDF$bg[statDF$bg == "RT_D3_Smooth.mm10.forModel.bedgraph"]                            <-  "Smooth"
  statDF$bg[statDF$bg == "RT_J185a_Myoblast_Int61896107_mm10.forModel.bedgraph"]           <-  "Myoblast"
  statDF$bg[statDF$bg == "RT_L1210_Lymphoblastoid_Ext73945012_mm10.forModel.bedgraph"]     <-  "LCL"

  ## set order (Mouse)
  bgOrder <- c("2C",
               "Sertoli",
               "4C",
               "MeiS(a)",
               "MeiS(b)",
               "MeiS(c)",
               "MeiS(SPO11)",
               "MeiS(wtLM)",
               "MeiS(S)a",
               "MeiS(S)b",
               "MeiS(H/D)",
               "MeiS(T1)",
               "MeiS(T2)a",
               "MeiS(T2)b",
               "MeiS(T3)",
               "MeiS(T4)",
               "MeiS(BxC)a",
               "MeiS(BxC)b",
               "MeiS(BxC)c",
               "SpgB(a)",
               "SpgB(b)",
               "SpgI(a)",
               "SpgI(b)",
               "SpgI(SPO11)",
               "SpgI(wtLM)",
               "SpgU(a)",
               "SpgU(b)",
               "SpgU(SPO11)",
               "SpgU(wtLM)",
               "Bcell",
               "CD8",
               "ESC(a)",
               "ESC(b)",
               "RPE1a",
               "RPE1b",
               "RPE1c",
               "BG02ESCa",
               "BG02ESCb",
               "Liver",
               "Myoblast",
               "Lymphoblastoid",
               "H7 ESC",
               "H9 NPC",
               "IMR90 Lung Fibro",
               "U2OS Bone Epi",
               "LCL",
               "Mammary",
               "MEF",
               "NPC",
               "PGC",
               "PRE-B",
               "Smooth",
               "SSC");

  ## Set ori order
  oriOrder <- c("OriSSDS(u)",
                "OriSSDS(hi)",
                'OriSSDS(WT1)',
                'OriSSDS(ESC)', 
                "OriSSDS(r)",
                "ATACSeq",
                "ATACSeq_at_CGI",
                "ATACSeq at CGI (Spg)",
                "ATACSeq no CGI (Spg)",
                "ATACSeq_no_CGI",
                "ATACSeq at CGI (LCL)",
                "ATACSeq no CGI (LCL)",
                "CGI",
                "CGI_no_ATACSeq",
                "CGI no ATACSeq (Spg)")

  ## Label "strength users"
  statDF$usingStr                    <- 'NO strength'
  statDF$usingStr[statDF$useStr]     <- 'Use strength'

  ## set order
  statDF$bg                         <- factor(statDF$bg,levels=bgOrder)
  statDF$ori                        <- factor(statDF$ori,levels=rev(oriOrder))

  ## Lable Meiotic samples
  statDF$cellType                                   <- 'Non-meiotic'
  statDF$cellType[grep(pattern = 'MeiS',statDF$bg)] <- 'Meiotic'

  ## Get the R2 cutoff for "best models"
  statAll <- aggregate(statDF$R2,
                       by=list(bg=statDF$bg),
                       FUN=function(x){quantile(x,bestModelCutoff,na.rm=TRUE)})

  statPer <- aggregate(statDF$R2,
                       by=list(bg=statDF$bg,ori=statDF$ori),
                       FUN=function(x){quantile(x,bestModelCutoff,na.rm=TRUE)})

  names(statPer)[3] <- 'R2cutoffPer'
  statPer$R2cutoffAll <- statAll$x

  statDF <- plyr:::join(statDF,statPer,by=c('bg','ori'))

  nGood <- which(statDF$R2 >= statDF$R2cutoffPer)
  nBad  <- sample(x       = which(statDF$R2 < statDF$R2cutoffPer),
                  size    = round(length(statDF$R2[statDF$R2<statDF$R2cutoffPer])/50),
                  replace = FALSE)

  dfAll <- statDF[c(nGood,nBad),]

  ## Get the origins that yield the best model
  ## Use these ONLY for downstream analyses
  bestOriSet                         <- unique(statDF$ori[statDF$R2 == max(statDF$R2)])[1]

  bestDF                             <-  statDF[statDF$ori == bestOriSet,]

  ## Downsample the crappy models (for computational speed)
  n                                  <- sample(x       = which(bestDF$R2 < bestDF$R2cutoffAll),
                                               size    = (sum(bestDF$R2 >= bestDF$R2cutoffAll)*10),
                                               replace = FALSE)

  bestDF                             <- bestDF[c(which(bestDF$R2 >= bestDF$R2cutoffAll),n),]

  ## Label "Best" models
  bestDF$topModels                                   <- 'Other models'
  bestDF$topModels[bestDF$R2 >= bestDF$R2cutoffAll]  <- 'Best models'

  bestDF$topModels                                   <- factor(bestDF$topModels,
                                                               levels=c('Best models','Other models'))

  return(list(df     = bestDF,
              dfAll  = dfAll,
              bgOrd  = bgOrder,
              oriOrd = oriOrder))
}

#### Some plots (NOT USED)
genModelPlots <- function(df,oSample,oOri){
  ## Plot 1: R2 plot
  gR2 <- ggplot(df,aes(fill=topModels,x=R2)) +
    geom_density(lwd=.4) +
    xlab(bquote('Model fit ('*R^2*")")) +
    scale_fill_manual(values=c('orange','grey50')) +
    theme(legend.position=c(0.5,1),
          legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key.size=unit(0.2,'cm')) +
    ylab('Model density')

  ## Plot 2: R2 by origins plot
  gOri <- ggplot(df,aes(x=topModels,y=oriPerMb,fill=topModels)) +
    geom_violin(adjust=2,trim = FALSE) + geom_boxplot(width=.1,fill='grey93') +
    scale_fill_manual(values=c('orange','grey50')) +
    theme(legend.position='none',
          legend.title=element_blank(),
          legend.justification=c(1,0),legend.background=element_blank(),
          legend.key.size=unit(0.2,'cm')) +
    xlab('') +
    ylab(bquote('Replisome density ('*Mb^-1*")")) +
    coord_flip()

  ## Plot 3: R2 by origins plot
  gTime <- ggplot(df,aes(x=topModels,y=repTime,fill=topModels)) +
    geom_violin(adjust=2,trim=FALSE) + geom_boxplot(width=.1,fill='grey93') +
    scale_fill_manual(values=c('orange','grey50')) +
    theme(legend.position='none',
          legend.title=element_blank(),
          legend.justification=c(1,0),legend.background=element_blank(),
          legend.key.size=unit(0.2,'cm')) +
    xlab('') +
    ylab(bquote('S-phase duration (cycles)')) +
    coord_flip(ylim=c(0,800))

  return(list(plotR2=gR2,
              plotOri=gOri,
              plotTime=gTime))
}

### IMPORT FCS FILE
importFCS <- function(fcsFile=NULL){
  fcm <- read.FCS(fcsFile,alter.names = TRUE)
  fcm <- as.data.frame((exprs(fcm)))
  return(fcm)
}

### GET SORTING DATA
getSort1Data <- function(fcsIn = '2019-09-11_WT_APC_PE_004.fcs',
                         limz = NULL) {
  dSort        <- importFCS(fcsIn)
  names(dSort)[names(dSort) == 'BV421.A'] <- 'DAPI'
  names(dSort)[names(dSort) == 'G575.A'] <- 'STRA8'
  names(dSort)[names(dSort) == 'B530.A'] <- 'SYCP3'
  names(dSort)[names(dSort) == 'R660.A'] <- 'DMRT1'

  dSort$DAPI2C  <- dSort$DAPI  > 80000 & dSort$DAPI < 115000
  dSort$DAPIRep  <- dSort$DAPI  > 115000 & dSort$DAPI < 145000
  dSort$DAPI4C  <- dSort$DAPI  > 145000 & dSort$DAPI < 220000
  dSort$DAPIok  <- dSort$DAPI  > 30000 & dSort$DAPI < 220000 & dSort$BV421.W < 88000
  dSort$STRA8ok  <- dSort$STRA8  > 1100
  dSort$DMRT1ok <- dSort$DMRT1 > 1000
  dSort$SYCP3ok <- dSort$SYCP3 > 800


  dSort$DSok    <- dSort$DAPIok

  dSort$allOK             <- 'Other'
  dSort$allOK[dSort$DSok] <- 'DAPI_Gated'
  dSort$ploidy <- '1C'
  dSort$ploidy[dSort$DAPI2C] <- '2C'
  dSort$ploidy[dSort$DAPIRep] <- 'Rep'
  dSort$ploidy[dSort$DAPI4C] <- '4C'

  dSortOK <- dSort[dSort$DSok,]

  return(list(all = dSort, ok=dSortOK))
}

drawCC_v_RR <- function(inDF,col1='SSDSaa1',col2='SSDSaa2',name='test'){
  ccDF <- data.frame(win=c(1,5,10,20,40,60,80,100,200,500,1000,1250,1500,1750,2000),ccOne=0,ccTwo=0)

  inDF$shuf1 <- sample(inDF[[col2]])
  inDF$shuf2 <- sample(inDF[[col2]])
  inDF$shuf3 <- sample(inDF[[col2]])
  inDF$shuf4 <- sample(inDF[[col2]])
  inDF$shuf5 <- sample(inDF[[col2]])
  inDF$shuf6 <- sample(inDF[[col2]])
  inDF$shuf7 <- sample(inDF[[col2]])
  inDF$shuf8 <- sample(inDF[[col2]])
  inDF$shuf9 <- sample(inDF[[col2]])
  inDF$shuf10 <- sample(inDF[[col2]])

  for (x in 1:length(ccDF$win)){
    wX <- ccDF$win[x]

    ccDF$cc12[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF[[col2]],wX,mean),method='spearman')
    ccDF$ccR1[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf1,wX,mean),method='spearman')
    # ccDF$ccR2[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf2,wX,mean))
    # ccDF$ccR3[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf3,wX,mean))
    # ccDF$ccR4[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf4,wX,mean))
    # ccDF$ccR5[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf5,wX,mean))
    # ccDF$ccR6[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf6,wX,mean))
    # ccDF$ccR7[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf7,wX,mean))
    # ccDF$ccR8[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf8,wX,mean))
    # ccDF$ccR9[x]  <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf9,wX,mean))
    # ccDF$ccR10[x] <- cor(rollapply(inDF[[col1]],wX,mean),rollapply(inDF$shuf10,wX,mean))
  }

  ccDF$name <- name;
  return(ccDF)
}


## Draw single-sim heatmap
plotSingleSimHeatmap <- function(modelData,
                                 csNum = 1,
                                 name = 'Meiotic',
                                 lowColor='green',
                                 midColor='black',
                                 highColor='violet'){
  
  #################################
  repPerTime <- data.frame(time=0:500,rep=0)
  for (i in 0:500){repPerTime$rep[(i+1)] <- sum(modelData$mod$repPerCell <= i)}
  
  repPerTime$pc <- repPerTime$rep / dim(modelData$mod$repPerCell)[1] /100 * 100
  
  ## Label 90% point as this is where we infer S-phase duration from
  pc90 <- which(abs(repPerTime$pc-90) == min(abs(repPerTime$pc-90)))[1]
  
  repPerTime$pc90 <- FALSE
  repPerTime$pc90[pc90] <- TRUE
  
  repPerTime$type <- name
  
  ######################## Heatmaps
  coords <- modelData$mod$model$whatCS == paste0('chr',csNum=csNum)
  
  mPerCell <- modelData$mod$repPerCell[coords,1:100]
  
  mPerCell[mPerCell == 0] <- NA
  mPerCell[mPerCell >repPerTime$time[pc90]] <- repPerTime$time[pc90]
  
  dfCells <- as.data.frame(mPerCell)
  
  dfCells$pos <- modelData$mod$model$trueCoord[coords]
  
  mdfCells <- reshape2:::melt.data.frame(dfCells,id.vars = 'pos',variable.name = 'cellid', value.name = 'RT')
  
  mdfCells$type <- name
  mdfCells$midPoint <- repPerTime$time[pc90]/2
  
  gHeatMapPerCell <- 
    ggplot(mdfCells) + geom_tile(aes(x=pos/1000000,
                                     y=cellid,
                                     fill=RT),
                                 na.rm=TRUE) + 
    scale_fill_gradient2(paste0(name,': simulated RT'),
                         low=lowColor,
                         mid=midColor,
                         high=highColor,
                         midpoint=mdfCells$midPoint[1],
                         na.value='white',
                         breaks=c(max(mPerCell,na.rm=TRUE)*.1,max(mPerCell,na.rm=TRUE)*.9),
                         labels=c('Early','Late'))+ 
    coord_cartesian(expand=FALSE) + 
    ylab('Modelled haploid genomes') + 
    xlab(paste0('Position on chromosome ',csNum,' (Mb)')) + 
    theme(legend.position='top',
          legend.direction = 'horizontal',
          legend.key.width=unit(2,'cm'),
          legend.key.height=unit(.3,'cm'),
          axis.text.y=element_blank(),
          legend.margin=margin(c(0,0,0,0)))
  
  return(list(data=mdfCells,
              fig=gHeatMapPerCell))
}

## Plot replicated DNA VS replication time
plotRepPerTime <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  
  repPerTime <- data.frame(time=0:500,rep=0)
  for (i in 0:500){repPerTime$rep[(i+1)] <- sum(modelData$mod$repPerCell <= i)}
  
  repPerTime$pc <- repPerTime$rep / dim(modelData$mod$repPerCell)[1] /100 * 100
  
  ## Label 90% point as this is where we infer S-phase duration from
  pc90 <- which(abs(repPerTime$pc-90) == min(abs(repPerTime$pc-90)))[1]
  
  repPerTime$pc90 <- FALSE
  repPerTime$pc90[pc90] <- TRUE
  
  repPerTime$type <- name
  
  gRepRateByTime <- ggplot(repPerTime,aes(x=time,
                                          y=pc)) + 
    geom_line(lwd=.3) + 
    geom_hline(yintercept=repPerTime$pc[which(repPerTime$pc90)],lwd=.2,color='magenta',lty='dashed') + 
    geom_vline(xintercept=repPerTime$time[which(repPerTime$pc90)],lwd=.2,color='magenta',lty='dashed')+ 
    ylab('Replicated (%)') + 
    xlab('Time (Simulation cycles)') + 
    annotate(geom='text',
             label=paste0('90% replicated'),
             size=8*5/14,
             x=10,y=repPerTime$pc[which(repPerTime$pc90)],
             hjust=0,vjust=1.4,
             color='magenta')
  
  if (returnFig){
    return(list(data=repPerTime,
                fig=gRepRateByTime))
  }else{
    return(list(data=repPerTime))
  }
}

## Plot active forks VS replication time
plotForksUsed <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  dfForks <- modelData$mod$activeForks$df %>% 
    group_by(time,cell) %>% 
    summarise(nf=sum(nForks),
              .groups='keep') %>% 
    group_by(time) %>% 
    summarise(nf=mean(nf)) %>%
    mutate(time = as.numeric(time)) %>% 
    dplyr:::filter(time <=500)
  
  gActiveForksByTime <- ggplot(dfForks,aes(x=as.numeric(time),y=nf)) + geom_point(alpha=1,size=.2) +
    ylab('Active forks (#)') + 
    xlab('Time (Model cycles)')
  
  dfForks$type <- name
  
  if (returnFig){
    return(list(data=dfForks,
                fig=gActiveForksByTime))
  }else{
    return(list(data=dfForks))
  }
}

## Plot origin firing frequency
plotOriFiringFrequency <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  
  ## GET DATA 
  dfCellOris <- data.frame(cell=1:100,nOris=0)
  
  for (n in 1:100){
    dfCellOris$nOris[n] <- dim(modelData$mod$forkDets[[paste0("cell",n)]])[1]/2
    if (n==1){
      oriID       <- modelData$mod$forkDets[[paste0("cell",n)]]$ori
    }else{
      oriID       <- c(oriID       ,modelData$mod$forkDets[[paste0("cell",n)]]$ori)
    }
  }
  
  ## ORI USAGE PLOT ############################################################
  dfOriFreq      <- as.data.frame(oriID) %>% group_by(oriID) %>% count() %>% mutate(freq=n/2) 
  dfOriUsage     <- dfOriFreq %>% group_by(freq) %>% count() %>% mutate(pc=n/sum(modelData$mod$model$ori>0)*100)
  
  dfOriUsage$medianFreq   <- median(dfOriFreq$freq)
  
  gOriUsage <- ggplot(dfOriUsage,aes(x=freq,y=pc)) + 
    geom_line(lwd=.3,color='grey60') + 
    geom_point(size=.3) + 
    ylab('Origins (%)') + 
    xlab('Firing frequency\n(Percent of simulations)') + 
    geom_vline(xintercept=dfOriUsage$medianFreq[1],
               lwd=.2,
               lty='dashed',
               color='magenta') +
    annotate(geom='text',
             label=paste0('median = ',dfOriUsage$medianFreq[1],' %'),
             size=8*5/14,
             x=dfOriUsage$medianFreq[1],y=Inf,
             hjust=-.1,vjust=1,
             color='magenta')
  
  dfOriUsage$type <- name
  
  if (returnFig){
    return(list(data=dfOriUsage,
                fig=gOriUsage))
  }else{
    return(list(data=dfOriUsage))
  }
}

## Plot origin used per sim
plotOrisPerSim <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  
  ## GET DATA 
  dfCellOris <- data.frame(cell=1:100,nOris=0)
  
  for (n in 1:100){
    dfCellOris$nOris[n] <- dim(modelData$mod$forkDets[[paste0("cell",n)]])[1]/2
    if (n==1){
      oriID       <- modelData$mod$forkDets[[paste0("cell",n)]]$ori
    }else{
      oriID       <- c(oriID ,modelData$mod$forkDets[[paste0("cell",n)]]$ori)
    }
  }
  
  ## ORI USAGE PLOT ############################################################
  gOrisUsed <- ggplot(dfCellOris,aes(x=nOris)) + 
    geom_histogram(bins = 5,fill='grey70',color='grey40') + 
    xlab('Origins fired (#)') + 
    ylab('Simulations (%)') + 
    annotate(geom='text',
             label=paste0('median:\n',format(round(median(dfCellOris$nOris)),big.mark = ",")," origins"),
             size=8*5/14,
             x=Inf,y=Inf,
             hjust=1,vjust=1,
             color='magenta')
  
  dfCellOris$type <- name
  
  if (returnFig){
    return(list(data=dfCellOris,
                fig=gOrisUsed))
  }else{
    return(list(data=dfCellOris))
  }
}

## Plot origin usage
plotOriUsage <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  
  ## GET DATA 
  dfCellOris <- data.frame(cell=1:100,nOris=0)
  
  for (n in 1:100){
    dfCellOris$nOris[n] <- dim(modelData$mod$forkDets[[paste0("cell",n)]])[1]/2
    if (n==1){
      vForkDist  <- modelData$mod$forkDets[[paste0("cell",n)]]$distTravelled
    }else{
      vForkDist  <- c(vForkDist  ,modelData$mod$forkDets[[paste0("cell",n)]]$distTravelled)
    }
  }
  
  ## FORK DISTANCE PLOT ########################################################
  numForks <- length(vForkDist)
  
  nStep          <- modelData$mod$params$repStep*1000
  medianForkDist <- median (vForkDist) * nStep / 1e6 ## Convert steps to Mb
  
  dfForkDist     <- as.data.frame(vForkDist) %>% 
    group_by(vForkDist) %>% 
    count() %>% 
    mutate(pc=n/numForks*100, 
           forkDist=vForkDist * nStep / 1e6) 
  
  dfForkDist$medianDist <- medianForkDist
  dfForkDist$type       <- name
  
  gForkDist      <- ggplot(dfForkDist,aes(x=forkDist,y=pc)) + 
    geom_line(lwd=.3,color='grey60') + 
    geom_point(size=.3) + 
    geom_vline(xintercept=dfForkDist$medianDist[1],lwd=.2,lty='dashed',color='magenta') + 
    geom_text(aes(x=medianDist,
                  label=paste0('\nMedian = ',medianDist,' Mb')),
              y=Inf,
              check_overlap=TRUE,
              size=8*5/14,
              color='magenta',
              hjust=-.2,vjust=1) + 
    ylab('Forks (%)') + 
    xlab('Distance replicated (Mb)')
  
  if (returnFig){
    return(list(data=dfForkDist,
                fig=gForkDist))
  }else{
    return(list(data=dfForkDist))
  }
}

## Plot origin usage
plotNewOriFiring <- function(modelData,name = 'Meiotic',returnFig=FALSE){
  
  ## GET DATA 
  dfCellOris <- data.frame(cell=1:100,nOris=0)
  
  for (n in 1:100){
    dfCellOris$nOris[n] <- dim(modelData$mod$forkDets[[paste0("cell",n)]])[1]/2
    if (n==1){
      vTimeFired <- modelData$mod$forkDets[[paste0("cell",n)]]$timeFired
    }else{
      vTimeFired <- c(vTimeFired ,modelData$mod$forkDets[[paste0("cell",n)]]$timeFired)
    }
  }
  
  ## NEW FIRING ORIGINS PLOT ###################################################
  numOrigins <- length(vTimeFired)/2 ## Num forks/2
  
  dfTimeFired     <- as.data.frame(vTimeFired) %>% 
    group_by(vTimeFired) %>% 
    count() %>% 
    mutate(n=n/2,
           percell=vTimeFired/100)
  
  gOriFiringRate <- ggplot(dfTimeFired,aes(x=vTimeFired,y = n/100)) + 
    scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),labels=c(0.01,0.1,1,10,100,1000)) + 
    geom_point(size=.3) + 
    geom_line(lwd=.3,color='grey60') + 
    geom_hline(yintercept=median(dfTimeFired$n/100),
               lty='dashed',
               lwd=.3,
               color='magenta') + 
    xlab('Time fired (cycle)') + 
    annotation_logticks(sides='l',
                        long  = unit(0.1,'cm'),
                        mid   = unit(0.05,'cm'),
                        short = unit(0.05,'cm'),
                        size=.2) +
    ylab('New origins fired\n(Mean count per simulation)') +
    annotate(geom='text',
             x=10,
             y=median(dfTimeFired$n/100),
             label=paste0('\n\nMedian = ',round(median(dfTimeFired$n/100),1),' ori per cycle'),
             size=8*5/14,
             color='magenta',
             hjust=0,vjust=1)
  
  dfTimeFired$type <- name
  
  if (returnFig){
    return(list(data=dfTimeFired,
                fig=gOriFiringRate))
  }else{
    return(list(data=dfTimeFired))
  }
}

## Plot RT vs model
#################################################
plotRTvModel <- function(dfMod,chrom2use="chr1",myName="tst",altTime=NULL){
  myModel <- dfMod$mod
  
  theme7point()
  
  if (chrom2use == 'all'){
    myRT                          <- esc$mod$experimentalRT$smooth
    # myModel$params$chrom          <- 'All'
    # myModel$model                 <- myModel$model[myModel$csPos,]
    # myName                        <- paste0(myName,'_All')
  }else{
    myModel$csPos                 <- myModel$model$whatCS == chrom2use
    myModel$params$chrom          <- chrom2use
    myModel$model                 <- myModel$model[myModel$csPos,]
    myName                        <- paste0(myName,'_',chrom2use)
    myRT                          <- myModel$experimentalRT$smooth[myModel$csPos]
  }
  
  if(is.null(altTime)){
    repTime <- dfMod$figData$bestStat$time
  }else{
    repTime <- altTime
  }

  mStat         <- testModelVData(myModel,
                                  myRT,
                                  dfMod$figData$bestStat$pcRep,
                                  repTime,
                                  noPlot=FALSE,
                                  normBoth=TRUE,
                                  smoothLevel=5,
                                  titleType = myName,
                                  upColor='#7fbf7b',
                                  downColor='#af8dc3')
  
  return (mStat$fig)
}

## Plot RT vs model
#################################################
plotRTvModel_wholeGenome <- function(dfMod){

  mStat         <- testModelVData(dfMod$mod,
                                  dfMod$mod$experimentalRT$smooth,
                                  dfMod$figData$bestStat$pcRep,
                                  dfMod$figData$bestStat$time,
                                  noPlot=FALSE,
                                  normBoth=TRUE,
                                  smoothLevel=5,
                                  titleType = 'NA',
                                  upColor='#7fbf7b',
                                  downColor='#af8dc3')
  
  return (mStat$fig)
}

## Make all model stats
makeAllStatsFigs <- function(modelFile,sampleName,chrom=1,modelData=NULL){
  if (!is.null(modelData)){
    rm('modelData'); gc()
    load(file = modelFile)
  }
  
  retList <- list()
  retList$lstHM         <- plotSingleSimHeatmap(modelData,chrom,sampleName)
  retList$lstRepRate    <- plotRepPerTime(modelData,sampleName)
  retList$lstForksUsed  <- plotForksUsed(modelData,sampleName)
  retList$lstOriFiring  <- plotOriFiringFrequency(modelData,sampleName)
  retList$lstOriPerSim  <- plotOrisPerSim(modelData,sampleName)
  retList$lstOriUsage   <- plotOriUsage(modelData,sampleName)
  retList$lstNewOriRate <- plotNewOriFiring(modelData,sampleName)
  
  # retList$gX6 <- ggarrange(retList$lstRepRate$fig,retList$lstOriPerSim$fig,
  #                          retList$lstForksUsed$fig,retList$lstOriFiring$fig,
  #                          retList$lstOriUsage$fig,retList$lstNewOriRate$fig,
  #                          ncol=2,nrow=3,
  #                          align='hv',
  #                          labels=c('B','D','C','E','F','G'),
  #                          font.label = list(size=8,face='bold'),
  #                          hjust=0,vjust=1)
  # 
  # retList$fig <- ggarrange(retList$lstHM$fig,retList$gX6,
  #                          ncol=1,heights=c(2,3),
  #                          labels=c('A',''),
  #                          font.label = list(size=8,face='bold'),
  #                          hjust=0,vjust=1)
  
  rm('modelData'); gc()
  return (retList)
}


