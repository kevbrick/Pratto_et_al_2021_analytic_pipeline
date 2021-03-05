## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')
source('accessoryFiles/scripts/R/simulateReplicationTiming_allCS.R')

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggridges)

# FUNCTIONS FOR FIGURE 4 #######################################################
plotRTvModel <- function(dfMod,chrom2use="chr1",myName="tst",altTime=NULL){
  myModel <- dfMod$mod
  
  theme7point()
  
  myModel$csPos                 <- myModel$model$whatCS == chrom2use
  myModel$params$chrom          <- chrom2use
  myModel$model                 <- myModel$model[myModel$csPos,]
  myName                        <- paste0(myName,'_',chrom2use)
  
  if(is.null(altTime)){
    repTime <- dfMod$figData$bestStat$time
  }else{
    repTime <- altTime
  }
  
  mStat         <- testModelVData(myModel,
                                  myModel$experimentalRT$smooth[myModel$csPos],
                                  dfMod$figData$bestStat$pcRep,
                                  repTime,
                                  noPlot=FALSE,
                                  normBoth=TRUE,
                                  smoothLevel=5,
                                  titleType = 'Internal',
                                  upColor='#7fbf7b',
                                  downColor='#af8dc3')
  
  return (mStat$fig)
}

renameSamplesAndSplit <- function(df){
  
  df$o <- 'Others'
  df$o[df$ori == 'OriSSDS(hi)']       <- 'Ori (hi-conf)'
  df$o[df$ori == 'OriSSDS(u)']        <- 'Ori (union)'
  df$o[df$ori == 'OriSSDS(r)']        <- 'Ori (random)'
  df$o[df$ori == 'OriSSDS(ESC)']      <- 'Ori (from WT 1 only)'
  df$o[df$ori == 'OriSSDS(WT1)']      <- 'Ori (from ESC only)'
  df$o[df$ori == 'ATACSeq_at_CGI']    <- 'AS+CGI'
  df$o[df$ori == 'ATACSeq_no_CGI']    <- 'AS-CGI'
  df$o[df$ori == 'CGI_no_ATACSeq']    <- 'CGI-AS'
  df$o[df$ori == 'Mechali Ori (ESC)'] <- 'ESC ori (Cayrou)'
  
  df <- df[df$o %in% c('Ori (hi-conf)',
                       'Ori (union)',
                       'Ori (random)',
                       'Ori (from WT 1 only)',
                       'Ori (from ESC only)',
                       'ESC ori (Cayrou)',
                       'AS+CGI',
                       'AS-CGI',
                       'CGI-AS'),]
  
  oriOrd <- c('Ori (hi-conf)',
              'Ori (union)',
              'Ori (from WT 1 only)',
              'Ori (from ESC only)',
              'ESC ori (Cayrou)',
              'AS+CGI',
              'AS-CGI',
              'CGI-AS',
              'Ori (random)')
  
  bgOrd <- c("Bcell",
             "LCL",
             "Myoblast",
             "CD8",
             "PGC",
             "SSC",
             "ESC(b)",
             "ESC(a)",
             "SpgU(SPO11)",
             "SpgU(wtLM)",
             "SpgU(b)",
             "SpgU(a)",
             "SpgI(SPO11)",
             "SpgI(wtLM)",
             "SpgI(b)",
             "SpgI(a)",
             "SpgB(b)",
             "SpgB(a)",
             "MeiS(SPO11)",
             "MeiS(wtLM)",
             "MeiS(c)",
             "MeiS(b)",
             "MeiS(a)",
             "2C",
             "4C",
             "MeiS(T1)",
             "MeiS(T2)a",
             "MeiS(T2)b",
             "MeiS(T3)",
             "MeiS(T4)")
  
  df$bg  <- factor(df$bg,levels=bgOrd)
  df$o   <- factor(df$o,levels=oriOrd)
  df$ori <- df$o
  
  df <- df[(df$bg %in% c('MeiS(a)','MeiS(b)','MeiS(c)','MeiS(wtLM)','MeiS(SPO11)',
                         'MeiS(T1)','MeiS(T2)a','MeiS(T2)b','MeiS(T3)','MeiS(T4)',
                         'SpgB(a)','SpgB(b)',
                         'SpgI(a)','SpgI(b)','SpgI(wtLM)','SpgI(SPO11)',
                         'SpgU(a)','SpgU(b)','SpgU(wtLM)','SpgU(SPO11)',
                         '2C','4C','Sertoli',
                         'ESC(a)','ESC(b)','SSC','PGC',
                         'CD8','Myoblast','LCL','Bcell')),]
  
  df <- df[!is.na(df$o),]
  df <- df[!is.na(df$bg),]
  
  ## For heatmaps
  dfPer <- df[df$R2>=df$R2cutoffPer,]
  
  dfPer <- dfPer[(dfPer$bg %in% c('2C','MeiS(a)','MeiS(b)','MeiS(c)',
                                  'SpgB(a)','SpgB(b)',
                                  'SpgI(a)','SpgI(b)',
                                  'SpgU(a)','SpgU(b)',
                                  'ESC(a)','ESC(b)','SSC','PGC',
                                  'CD8','Myoblast','LCL','Bcell')),]
  
  return(list(dfA   = df,
              dfPer = dfPer))
}

drawOri4ModelHeatmap <- function(dfPer,loCol='dodgerblue1',hiCol='red'){
  
  pc100 <- function(x){
    if (x >= 0.1){
      return(x*100)
    }else{
      return(paste0("0",x*100))
    }
  }
  dfPlot <- dfPer %>% mutate(ori=o, sample=bg) %>%
    group_by(sample,ori) %>%
    summarise(meanR2      = round(mean(R2),2),
              label       = paste0("0.",pc100(meanR2)),
              labelNoZero = paste0(".",pc100(meanR2))) %>%
    ungroup() %>% mutate(szR2=standardize0to1(0.5*(1-meanR2))) %>%
    select(sample,ori,meanR2,label,labelNoZero,szR2)
  
  gSized <- ggplot(dfPlot) + 
    geom_tile(aes(x=ori,y=sample,fill=meanR2,size=1.4*szR2),color='white') +
    scale_fill_gradient2(bquote("Model fit (Mean "*R^2*")"),
                         low=loCol,
                         mid='grey90',
                         high=hiCol,
                         midpoint=0.3) +
    xlab('') +
    ylab('') +
    theme(legend.position='top',
          legend.direction = 'horizontal',
          legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
          legend.key.size=unit(0.1,'cm'),
          legend.key.width=unit(0.3,'cm')) + 
    guides(size = FALSE)
  
  gLblNoZero <- ggplot(dfPlot) + 
    geom_tile(aes(x=ori,y=sample,fill=meanR2),color='white') +
    geom_text(aes(x=ori,y=sample,label=labelNoZero),size=7*5/14) +
    scale_fill_gradient2(bquote("Model fit (Mean "*R^2*")"),
                         low=loCol,
                         mid='grey90',
                         high=hiCol,
                         midpoint=0.3) +
    xlab('') +
    ylab('') +
    theme(legend.position='top',
          legend.direction = 'horizontal',
          legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
          legend.key.size=unit(0.1,'cm'),
          legend.key.width=unit(0.3,'cm')) + 
    guides(size = FALSE)
  
  gLbl <- ggplot(dfPlot) + 
    geom_tile(aes(x=ori,y=sample,fill=meanR2),color='white') +
    geom_text(aes(x=ori,y=sample,label=label),size=7*5/14) +
    scale_fill_gradient2(bquote("Model fit (Mean "*R^2*")"),
                         low=loCol,
                         mid='grey90',
                         high=hiCol,
                         midpoint=0.3) +
    xlab('') +
    ylab('') +
    theme(legend.position='top',
          legend.direction = 'horizontal',
          legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
          legend.key.size=unit(0.1,'cm'),
          legend.key.width=unit(0.3,'cm')) + 
    guides(size = FALSE)
  
  gNoLbl <- ggplot(dfPlot) + 
    geom_tile(aes(x=ori,y=sample,fill=meanR2),color='white',lwd=.3) +
    scale_fill_gradient2(bquote("Model fit (Mean "*R^2*")"),
                         low=loCol,
                         mid='grey90',
                         high=hiCol,
                         midpoint=0.3) +
    xlab('') +
    ylab('') +
    theme(legend.position='top',
          legend.direction = 'horizontal',
          legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
          legend.key.size=unit(0.1,'cm'),
          legend.key.width=unit(0.3,'cm')) + 
    guides(size = FALSE)
  
  commonTheme <- theme(plot.margin=margin(0,0,0,0,'cm'),
                       legend.text=element_text(size=7),
                       axis.text.x=element_text(angle=45,hjust=1)) 
  
  return(list(gSz=gSized + commonTheme,
              gNZ=gLblNoZero + commonTheme,
              gLbl=gLbl + commonTheme,
              gNoLbl=gNoLbl + commonTheme))
}

drawSpo11KOSortingSchema <- function(pointSz = 0.2, pointA  = 0.2){
  library(data.table)
  
  theme7point()
  
  fcsIn ='Specimen_001_Spo-11_001.fcs'
  limz = NULL
  
  dSort        <- importFCS(fcsIn)
  names(dSort)[names(dSort) == 'BV421.A'] <- 'DAPI'
  names(dSort)[names(dSort) == 'G575.A'] <- 'STRA8'
  names(dSort)[names(dSort) == 'B530.A'] <- 'SYCP3'
  names(dSort)[names(dSort) == 'R660.A'] <- 'DMRT1'
  
  dSort$DAPI2C  <- dSort$DAPI  > 60000 & dSort$DAPI < 95000
  dSort$DAPIRep <- dSort$DAPI  > 95000 & dSort$DAPI < 125000
  dSort$DAPI4C  <- dSort$DAPI  > 125000 & dSort$DAPI < 200000
  dSort$DAPIok  <- dSort$DAPI  > 40000 & dSort$DAPI < 180000 & dSort$BV421.W < 110000
  
  dSort$DSok    <- dSort$DAPIok
  
  dSort$allOK                 <- 'Unclassified'
  dSort$allOK[dSort$DSok]     <- 'DAPI_Gated'
  dSort$ploidy                <- '1C'
  dSort$ploidy[dSort$DAPI2C]  <- '2C'
  dSort$ploidy[dSort$DAPIRep] <- 'Rep'
  dSort$ploidy[dSort$DAPI4C]  <- '4C'
  
  dSortOK <- dSort[dSort$DSok,]
  
  
  # colz  <- c("SpgB" = "#009E73", 
  #            "SpgU"="#CC79A7", 
  #            "SpgI"="#D55E00", 
  #            "MeiS" = "#444444",
  #            "Unclassified" = "#cccccc", 
  #            "Other1" = "#E69F00", 
  #            "Other2" = "#56B4E9", 
  #            "Other3" = "#0072B2",
  #            "Other4" = "#F0E442")
  # 
  colz  <- c("SpgB" = "#D58C00", 
             "SpgU"="#D57000", 
             "SpgI"="#D55E00", 
             "MeiS" = "#CC79A7",
             "Unclassified" = "#cccccc", 
             "Other1" = "#E69F00", 
             "Other2" = "#56B4E9", 
             "Other3" = "#0072B2",
             "Other4" = "#F0E442")
  
  dRep <- dSortOK[dSortOK$DAPIRep,]
  d4C  <- dSortOK[dSortOK$DAPI4C,]
  d2C  <- dSortOK[dSortOK$DAPI2C,]
  nTot <- length(dSortOK$DAPI)
  nW   <- 750
  
  dDAPI    <- as.data.frame(table(round(dSortOK$DAPI/nW)))
  dDAPIRep <- as.data.frame(table(round(dSortOK$DAPI[dSortOK$DAPIRep]/nW)))
  
  dDAPI$DAPI    <- as.numeric(as.character(dDAPI$Var1))*nW/10000
  dDAPIRep$DAPI <- as.numeric(as.character(dDAPIRep$Var1))*nW/10000
  
  dDAPI$pc      <- dDAPI$Freq/nTot*100
  dDAPIRep$pc   <- dDAPIRep$Freq/nTot*100
  
  dDAPI         <- dDAPI[dDAPI$DAPI < 19.5,]
  dDAPIRep      <- dDAPIRep[dDAPIRep$DAPI < 19.5,]
  
  dlDAPI        <- dDAPI[dDAPI$DAPI > 7.5 & dDAPI$DAPI < 19.5,]
  dlDAPIRep     <- dDAPIRep[dDAPIRep$DAPI > 7.5 & dDAPIRep$DAPI < 19.5,]
  
  gALL <- ggplot() + geom_area(data=dDAPI,aes(x=DAPI,y=pc),color='grey40',lwd=.2) +
    geom_area(data=dDAPIRep,aes(x=DAPI,y=pc),fill='red',color='firebrick',lwd=.2) +
    coord_cartesian(expand=FALSE) +
    xlab(bquote('DAPI intensity (x'*10^4*')')) +
    ylab('Nuclei (%)')
  
  dLabels <- data.frame(x=c(8,11,16),
                        y=c(0.05,1.0,0.05),
                        lbl=c('2C','2-4C','4C'))
  
  gDNAContent <- ggplot() +
    annotate(geom='rect',xmin=min(dDAPIRep$DAPI),xmax=max(dDAPIRep$DAPI),ymin=0.001,ymax=10,
             fill='firebrick',alpha=.2,color='firebrick',lwd=.3) +
    coord_cartesian(ylim=c(0.02,8),expand=FALSE) +
    geom_smooth(data=dDAPI,aes(x=DAPI,y=pc),color='grey40',lwd=.4,span=.2,se=FALSE) +
    xlab(bquote('DAPI intensity (x'*10^4*')')) +
    ylab('Nuclei (%)') + scale_y_log10(breaks=c(0.1,1)) +
    geom_text(data=dLabels,aes(x=x,y=y,label=lbl),size=7*5/14,fontface='bold') +
    annotation_logticks(sides='l',
                        long=unit(0.2,'cm'),
                        mid=unit(0.1,'cm'),
                        short=unit(0.1,'cm'),
                        size=.2)
  
  
  dRep$stage <- "Unclassified"
  dRep$stage[dRep$DMRT1 > 2000 & dRep$DMRT1 < 5000 & dRep$STRA8 > 600  & dRep$STRA8 < 3000 ] <- 'SpgU'
  dRep$stage[dRep$DMRT1 > 5000 & dRep$DMRT1 < 30000 & dRep$STRA8 < 20000  & dRep$STRA8 > 3100] <- 'SpgB'
  dRep$stage[dRep$DMRT1 < 3000 & dRep$DMRT1 > 200 & dRep$STRA8 < 20000  & dRep$STRA8 > 5000] <- 'MeiS'
  dRep$stage[dRep$DMRT1 > 5000 & dRep$DMRT1 < 30000 & dRep$STRA8 < 3000  & dRep$STRA8 >600 ] <- 'SpgI'
  #
  # fortxt <- data.frame(stage = c('Late MeiS/Leptotene','Spg', 'MeiS'),
  #                      H1 = c(1000, 25 , 25),
  #                      DMRT1 = c(12, 50000, 5000))
  dRep$stage <- factor(dRep$stage, levels=c('MeiS','SpgU','SpgB','SpgI','Unclassified'))
  
  
  
  fortxt <- data.frame(stage = c('MeiS','SpgU','SpgB','SpgI'),
                       STRA8 = c(15000, 1000, 15000, 600),
                       DMRT1 = c(2000, 3500, 20000, 20000))
  
  gSort <- ggplot(dRep,aes(x=STRA8,y=DMRT1,
                           fill=as.factor(stage),
                           color=as.factor(stage))) +
    geom_density2d(lwd=.05,alpha=0.8) +
    geom_point(shape=20,size=pointSz,alpha=pointA) +
    scale_x_log10(breaks=c(100,1000,10000), labels=fancy_scientific) +
    scale_y_log10(breaks=c(100,1000,10000), labels=fancy_scientific) +
    annotate(geom='rect',
             xmin=600,xmax = 3000,
             ymax=5000,ymin=2000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'SpgU'])+
    annotate(geom='rect',
             xmin=3100,xmax=20000,
             ymax=30000,ymin=5000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'SpgB'])+
    annotate(geom='rect',
             xmax=20000,xmin=5000,
             ymax=3000,ymin=200,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'MeiS'])+
    annotate(geom='rect',
             xmax=3000,xmin=600,
             ymax=30000,ymin=5000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'SpgI'])+
    annotation_logticks(sides='bl',
                        size=0.2,
                        short=unit(0.050,'cm'),
                        mid=unit(0.075,'cm'),
                        long=unit(0.100,'cm'))+
    scale_color_manual(values=colz) +
    theme(legend.title=element_blank(),
          legend.position='none',
          legend.justification=c(1,0),
          legend.direction = 'horizontal',
          legend.key.size = unit(.2,'cm'),
          legend.background=element_blank()) +
    xlab('STRA8 fluorescence') +
    ylab('DMRT1 fluorescence') +
    coord_cartesian(ylim=c(20,40000),xlim=c(400,30000))+
    annotate(geom='text',
             x=20000,y=200,
             label='MeiS',
             hjust=1,vjust=1.5,
             size=8*5/14,
             color=colz[names(colz) == 'MeiS']) +
    annotate(geom='text',
             x=600,y=2000,
             label='SpgU',
             hjust=0,vjust=1.5,
             size=8*5/14,
             color=colz[names(colz) == 'SpgU'])+
    annotate(geom='text',x=20000,y=30000,
             label='SpgB',
             hjust=1,vjust=-.5,
             size=8*5/14,
             color=colz[names(colz) == 'SpgB']) +
    annotate(geom='text',
             x=600,y=30000,
             label='SpgI',
             hjust=0,vjust=-0.5,
             size=8*5/14,
             color=colz[names(colz) == 'SpgI'])
  
  gSpo11Sort <- ggarrange(gDNAContent, gSort, labels = c("A","B"), font.label = list(size=8))
  
  return(gSpo11Sort)
}

plotRidges <- function(dR,pType,myAlpha=0.9){
  if (pType == 'main'){
    dR <- dR[(dR$bg %in% c('MeiS(a)','MeiS(b)','MeiS(c)',
                           'SpgB(a)','SpgB(b)','SpgI(a)','SpgI(b)','SpgU(a)','SpgU(b)',
                           'ESC(a)','ESC(b)','SSC','PGC',
                           'CD8','Myoblast','LCL','Bcell')),]
  }
  
  if (pType == 'suppSpo11'){
    dR <- dR[(dR$bg %in% c('MeiS(SPO11)','MeiS(wtLM)','MeiS(a)','MeiS(b)','MeiS(c)',
                           'SpgI(SPO11)','SpgI(wtLM)','SpgI(a)','SpgI(b)',
                           'SpgU(SPO11)','SpgU(wtLM)','SpgU(a)','SpgU(b)')),]
    
    #dR <- dR[dR$ori == 'Ori (hi-conf)',]
  }
  
  if (pType == 'all'){
    #dR <- dR[dR$ori == 'Ori (hi-conf)' & (dR$bg %in% c('MeiS(SPO11)','MeiS(wtLM)','MeiS(a)','MeiS(b)','MeiS(c)',
    dR <- dR[(dR$bg %in% c('MeiS(SPO11)','MeiS(wtLM)','MeiS(a)','MeiS(b)','MeiS(c)',
                           'SpgB(a)','SpgB(b)',
                           'SpgI(SPO11)','SpgI(wtLM)','SpgI(a)','SpgI(b)',
                           'SpgU(SPO11)','SpgU(wtLM)','SpgU(a)','SpgU(b)',
                           'ESC(a)','ESC(b)','SSC','PGC',
                           'CD8','Myoblast','LCL','Bcell')),]
  }
  
  palette5 <- c('#beaed4',
                '#de77ae',
                '#fdc086',
                '#ffff99',
                '#386cb0')
  
  palette4 <- c('#de77ae',
                'orange',
                '#555555',
                'yellow',
                '#bbbbbb')
  
  dR <- dR[dR$R2 >= dR$R2cutoffAll,]
  
  gR2 <- ggplot(dR,aes(y=bg,x=R2,fill=type)) +
    geom_density_ridges(alpha=myAlpha,bandwidth=0.022,lwd=.1) +
    theme(legend.position='none') + ylab('') +
    xlab(bquote("Model fit ("*R^2*')')) +
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values=palette4)
  
  dR <- dR[!(dR$bg %in% c('2C','4C','Sertoli')),]
  
  gOri <- ggplot(dR,aes(y=bg,x=oriPerMb*2,fill=type)) +
    geom_density_ridges(alpha=myAlpha,bandwidth=0.022,lwd=.1) +
    theme(legend.position='none') + ylab('') + xlab(bquote("Firing rate (Replisomes "*Mb^-1*')')) +
    coord_cartesian(expand=FALSE) +
    scale_fill_manual(values=palette4)
  
  mO  <- as.data.frame(aggregate(dR$oriPerMb,list(bg = dR$bg, type=dR$type),FUN=mean))
  mdO <- as.data.frame(aggregate(dR$oriPerMb,list(bg = dR$bg, type=dR$type),FUN=median))
  sdO <- as.data.frame(aggregate(dR$oriPerMb,list(bg = dR$bg, type=dR$type),FUN=sd))
  
  mR  <- as.data.frame(aggregate(dR$repTime,list(bg = dR$bg, type=dR$type),FUN=mean))
  mdR <- as.data.frame(aggregate(dR$repTime,list(bg = dR$bg, type=dR$type),FUN=median))
  sdR <- as.data.frame(aggregate(dR$repTime,list(bg = dR$bg, type=dR$type),FUN=sd))
  ssRT <- cbind(mdO,mO$x,sdO$x,mdR$x,mR$x,sdR$x)
  names(ssRT) <- c('bg','type','medianOriPerMb','meanOriPerMb','sdOriPerMb','medianRT','meanRT','sdRT')
  
  gRT <- ggplot(dR,aes(y=bg,x=repTime,fill=type)) +
    geom_density_ridges(alpha=myAlpha,bandwidth=22,lwd=.1) +
    theme(legend.position='none') +
    coord_cartesian(expand=FALSE) +
    xlab(bquote("S-phase duration (cycles)") ) +
    ylab('') +
    scale_fill_manual(values=palette4)
  
  return(list(gR2 = gR2,
              gOri = gOri,
              gRT = gRT,
              ssData=ssRT))
}

drawModelMetrics <- function(df,myAlpha=0.9){
  ### Draw VE->L plot first
  d <- df
  
  #MeiS(T2)b
  
  d <- d[(d$bg %in% c('MeiS(T1)',
                      'MeiS(T2)a',
                      'MeiS(T3)',
                      'MeiS(T4)')),]
  
  d$nm <- 'T1'
  d$nm[d$bg == 'MeiS(T2)a'] <- 'T2'
  d$nm[d$bg == 'MeiS(T2)b'] <- 'T2b'
  d$nm[d$bg == 'MeiS(T3)']  <- 'T3'
  d$nm[d$bg == 'MeiS(T4)']  <- 'T4'
  
  d$bg <- d$nm
  
  bgOrd <- rev(c('T1',
                 'T2',
                 'T2a',
                 'T3',
                 'T4'))
  
  d$bg <- factor(d$bg,levels=bgOrd)
  
  palette5 <- c('#beaed4',
                '#de77ae',
                '#fdc086',
                '#ffff99',
                '#386cb0')
  
  palette4 <- c('#de77ae',
                'orange',
                '#555555',
                'yellow',
                '#bbbbbb')
  
  d <- d[d$R2 >= d$R2cutoffAll,]
  
  gT14 <- ggplot(d,aes(y=bg,x=time,fill=type)) +
    geom_density_ridges(alpha=myAlpha,bandwidth=22,lwd=.1) +
    theme(legend.position='none') +
    coord_cartesian(expand=FALSE) +
    scale_x_continuous(breaks=c(0,75,150,225,300,375),labels=c(0,'',150,'',300,''))+
    xlab("Optimal run-time\n(cycles)") +
    ylab('') +
    scale_fill_manual(values=palette4)
  
  gSlen14 <- ggplot(d,aes(y=bg,x=repTime,fill=type)) +
    geom_density_ridges(alpha=myAlpha,bandwidth=22,lwd=.1) +
    theme(legend.position='none') +
    coord_cartesian(expand=FALSE) +
    xlab("S-phase duration\n(cycles)") +
    scale_x_continuous(breaks=c(200,300,400,500,600),labels=c(200,'',400,'',600))+
    ylab('') +
    scale_fill_manual(values=palette4)
  
  ## All others
  d <- df
  lrMain  <- plotRidges(d,'main')
  lrSpo11 <- plotRidges(d,'suppSpo11')
  lrAll   <- plotRidges(d,'all')
  
  ssX2   <- lrAll$ssData
  
  # rtMeiSWT <- mean(ssX2$medianRT[grep('MeiS\\([abc]',ssX2$bg)])
  # 
  # ssX2$rtVMeiSwt      <- round(ssX2$medianRT/rtMeiSWT*100)

  rtMeiSWT <- mean(ssX2$meanRT[grep('MeiS\\([abc]',ssX2$bg)])
  
  ssX2$rtVMeiSwt      <- round(ssX2$meanRT/rtMeiSWT*100)
  ssX2$medianOriPerMb <- round(ssX2$medianOriPerMb*2,2)
  ssX2$meanOriPerMb   <- round(ssX2$meanOriPerMb*2,2)
  ssX2$sdOriPerMb     <- round(ssX2$sdOriPerMb*2,2)
  ssX2$medianRT       <- round(ssX2$medianRT)
  ssX2$meanRT         <- round(ssX2$meanRT)
  ssX2$sdRT           <- round(ssX2$sdRT)
  
  gBars <- ggplot(ssX2,aes(x=bg,y=rtVMeiSwt,fill=type)) +
    geom_bar(stat='identity',lwd=.2,color='grey30') +
    coord_flip() +
    geom_hline(yintercept=100,lty='dashed',lwd=.3,color='red') +
    ylab('S-phase duration\n(relative to MeiS (%))') +
    xlab('') +
    scale_fill_manual(values=palette4) +
    theme(legend.key.size = unit('0.3','cm'),
          legend.title=element_blank())
  
  #---------------------- TABLE
  dfTable    <- ssX2 %>% select(-type, -rtVMeiSwt)
  
  gStatTab   <- ggtexttable(dfTable,
                            rows=NULL,
                            cols = c('Sample',
                                     'Replisomes\nper Mb\n(median)',
                                     'Replisomes\nper Mb\n(mean)',
                                     'Replisomes\nper Mb\n(s.d.)',
                                     'S-phase\nduration\n(median)',
                                     'S-phase\nduration\n(mean)',
                                     'S-phase\nduration\n(s.d.)'),
                            theme=ttheme(colnames.style = colnames_style(size = 6),
                                         rownames.style = rownames_style(size = 6),
                                         tbody.style = tbody_style(size = 7),
                                         padding = unit(c(.6, 1.2), "mm"))
  )
  
  return(list(gR2     = lrMain$gR2,
              gOri    = lrMain$gOri,
              gRT     = lrMain$gRT,
              ssMain  = lrMain$ssData,
              gSR2    = lrSpo11$gR2,
              gSOri   = lrSpo11$gOri,
              gSRT    = lrSpo11$gRT,
              ssSpo11 = lrSpo11$ssData,
              ssALL   = ssX2,
              gT14    = gT14,
              gT14S   = gSlen14,
              gRTrel  = gBars,
              dfTab   = dfTable,
              gTable  = gStatTab))
}

drawSortSections <- function(){
  sd1 <- getSort1Data()
  
  dOK <- sd1$ok
  
  dd <- aggregate(dOK$DAPI,by=list(dapi=round(dOK$DAPI/50)*50),FUN=function(x){sum(x>0)})
  
  dfBoxes <- data.frame(xlo=c(102500,120000,130000,140000)/10000,
                        xhi=c(110000,127500,137500,147500)/10000)
  
  dfText <- data.frame(x=c(4.85,9.2,(dfBoxes$xlo+dfBoxes$xhi)/2,17.5),
                       y=c(100,10,100,20,100,20,15),
                       lbl=c('1C','2C','T1','T2','T3','T4','4C'))
  
  gSort <- ggplot() + scale_y_log10() +
    geom_rect(data=dfBoxes,
              aes(xmin=xlo,xmax=xhi),
              ymin=-Inf,ymax=Inf,
              fill=alpha('#de77ae',0.4)) +
    geom_smooth(data=dd,
                aes(x=dapi/10000,y=x),
                span=.15,
                method='loess',
                color='grey50',
                lwd=.4) +
    geom_text(data=dfText,
              aes(x=x,y=y,label=lbl),
              size=7*5/14,fontface='bold') +
    xlab(bquote('DAPI (x'*10^5*")")) +
    ylab('Nuclei') +
    annotation_logticks(sides='l',
                        long = unit(0.2,'cm'),
                        mid =  unit(0.1,'cm'),
                        short = unit(0.1,'cm'),
                        size = 0.2)
  
  return(gSort)
}

drawModelFitTable <- function(){
  
  files <- list.files(path=".", pattern="*mm10.modelMetrics.Rdata$", full.names=TRUE, recursive=FALSE)
  
  for(i in 1:length(files)){
    load(files[i])
    
    if (exists('metricsDF')){
      metricsDF$df    <- rbind(metricsDF$df,    dfMetrics$df)
      metricsDF$dfAll <- rbind(metricsDF$dfAll, dfMetrics$dfAll)
    }else{
      metricsDF <- dfMetrics
    }
  }
  
  lstData <- renameSamplesAndSplit(metricsDF$dfAll)
  
  lstData$dfPer <- lstData$dfPer %>% dplyr:::filter(!grepl('Cayrou',ori))

  ## Draw heatmap table
  lstHMs <- drawOri4ModelHeatmap(lstData$dfPer)
  gHM <- lstHMs$gLabelled 
  
  return(list(fig=lstHMs$gLbl,
              figs=lstHMs,
              data=metricsDF$df))
}

drawRTEMLFig <- function(){
  
  dfEML <- fread('rep_v_rec_MM10.tab')
  
  cs2use <- 'chr13'
  
  dfPlotEML <- dfEML %>% mutate(T1=expRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10,
                                T2=expRT_MeiS_EARLY_S3_2to4C_SCP3STRA8_NA_mm10,
                                T3=expRT_MeiS_MID_S2_2to4C_SCP3yH2AX_NA_mm10,
                                T4=expRT_MeiS_LATE_S2_2to4C_SCP3yH2AX_NA_mm10) %>%
    dplyr:::filter(cs == cs2use) %>%
    select(from,T1,T2,T3,T4) %>% 
    reshape2:::melt.data.frame(id.vars = 'from',
                               variable.name = 'sample', 
                               value.name = 'RT') %>%
    mutate(RTneg=ifelse(RT>0,0,RT),
           RTpos=ifelse(RT<0,0,RT))
  
  gFig <- ggplot(dfPlotEML,aes(x=from/1000000,y=RT)) +
    geom_hline(yintercept=0,lwd=.3)+
    geom_ribbon(aes(ymin=0,ymax=RTpos),fill=alpha('darkolivegreen',.5))+
    geom_ribbon(aes(ymin=0,ymax=RTneg),fill=alpha('purple',.5))  +
    facet_grid(sample~.) +
    xlab(paste0('Position on chr ',cs2use,' (Mb)')) +
    ylab(bquote('Relative RT ; '*log[2]*'(RT/median)')) +
    scale_y_continuous(breaks=c(-0.4,0,0.4)) +
    theme(panel.border = element_rect(size=.5,color='grey30',fill=NA))
  
  return(list(fig=gFig,
              data=dfPlotEML))
  
}

# PLOT FIGURE ##################################################################
theme7point()

chrom <- 8

load(file = 'ESC_ALL_S1_2to4C_NA_NA_mm10.Rdata')
esc <- modelData
rm('modelData'); gc()

## Set heatmap themes
thmHeatmapESC <-  theme(legend.title=element_blank(),
                        legend.margin=margin(c(0,0,-0.3,0),unit='cm'),
                        legend.position=element_blank(),
                        axis.ticks.y = element_blank(),
                        plot.margin=margin(c(0,0,0,0),unit='cm'),
                        axis.text.x=element_blank())  

sampleName <- 'ESC'
esc$lstHM         <- plotSingleSimHeatmap(esc,chrom,sampleName)
esc$lstRepRate    <- plotRepPerTime(esc,sampleName)
esc$lstForksUsed  <- plotForksUsed(esc,sampleName)
esc$lstOriFiring  <- plotOriFiringFrequency(esc,sampleName)
esc$lstOriPerSim  <- plotOrisPerSim(esc,sampleName)
esc$lstOriUsage   <- plotOriUsage(esc,sampleName)
esc$lstNewOriRate <- plotNewOriFiring(esc,sampleName)
esc$meanRT        <- plotRTvModel(esc,"chr8","cs1ESC")
esc$genomeRT      <- plotRTvModel_wholeGenome(esc)

esc$gMeanSC <- grid.arrange(left=textGrob(label='ES cell S-phase',gp = gpar(fontface='bold',fontsize=7),hjust=0.5,vjust=1,rot=90),
                            ggarrange(esc$lstHM$fig + thmHeatmapESC + theme(legend.position='none')+
                                        ylab('scRT-Sim') + 
                                        xlab('') + 
                                        theme(axis.text.x=element_blank()),
                                      esc$meanRT + 
                                        coord_cartesian(ylim=c(-2.5,2.5),xlim=c(0,129.4),expand=FALSE)+
                                        ylab('normalized RT') + 
                                        scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(-2,'',0,'',2)) +
                                        theme(plot.margin=unit(c(0,0,0,0.3),'cm')),
                                      ncol=1,nrow=2,
                                      align='v',
                                      heights=c(1,1)))

esc$gMeanSCemptyHM <- grid.arrange(left=textGrob(label='ES cell S-phase',gp = gpar(fontface='bold',fontsize=7),hjust=0.5,vjust=1,rot=90),
                            ggarrange(ggplot() + thmHeatmapESC + theme(legend.position='none')+
                                        ylab('scRT-Sim') + 
                                        xlab('') + 
                                        theme(axis.text.x=element_blank()),
                                      esc$meanRT + 
                                        coord_cartesian(ylim=c(-2.5,2.5),xlim=c(0,129.4),expand=FALSE)+
                                        ylab('normalized RT') + 
                                        scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(-2,'',0,'',2)) +
                                        theme(plot.margin=unit(c(0,0,0,0.3),'cm')),
                                      ncol=1,nrow=2,
                                      align='v',
                                      heights=c(1,1)))

load(file = 'MeiS_ALL_S3_2to4C_STRA8_DMRT1_mm10.Rdata')
meiS <- modelData
rm('modelData'); gc()

sampleName <- 'MeiS'
meiS$lstHM         <- plotSingleSimHeatmap(meiS,chrom,sampleName)
meiS$lstRepRate    <- plotRepPerTime(meiS,sampleName)
meiS$lstForksUsed  <- plotForksUsed(meiS,sampleName)
meiS$lstOriFiring  <- plotOriFiringFrequency(meiS,sampleName)
meiS$lstOriPerSim  <- plotOrisPerSim(meiS,sampleName)
meiS$lstOriUsage   <- plotOriUsage(meiS,sampleName)
meiS$lstNewOriRate <- plotNewOriFiring(meiS,sampleName)
meiS$meanRT        <- plotRTvModel(meiS,"chr8","cs1meiS")
meiS$genomeRT      <- plotRTvModel_wholeGenome(meiS)

thmHeatmapMeiS <-  theme(legend.title=element_blank(),
                         legend.margin=margin(c(0,0,-0.3,0),unit='cm'),
                         legend.text=element_text(size=7,vjust=7),
                         legend.text.align = c(4,-3),
                         legend.key.width=unit(1.2,'cm'),
                         legend.key.height=unit(0.2,'cm'),
                         axis.ticks.y = element_blank(),
                         plot.margin=margin(c(-0.2,0,-0.2,0),unit='cm'),
                         axis.text.x=element_blank())   

meiS$gMeanSC <- grid.arrange(left=textGrob(label='Meiotic S-phase',gp = gpar(fontface='bold',fontsize=7),rot=90),
                             ggarrange(meiS$lstHM$fig + thmHeatmapMeiS + 
                                         theme(legend.position='none') +
                                         ylab('scRT-Sim') + 
                                         xlab('') + 
                                         theme(axis.text.x=element_blank()),
                                       meiS$meanRT + 
                                         coord_cartesian(ylim=c(-2.5,2.5),xlim=c(0,129.4),expand=FALSE) +
                                         xlab('') + 
                                         ylab('normalized RT') + 
                                         scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(-2,'',0,'',2)) +
                                         theme(plot.margin=unit(c(0,0,0.1,0.3),'cm'),
                                               axis.text.x=element_blank()),
                                       ncol=1,nrow=2,
                                       align='v',
                                       heights=c(1,1)))

meiS$gMeanSCemptyHM <- grid.arrange(left=textGrob(label='Meiotic S-phase',gp = gpar(fontface='bold',fontsize=7),rot=90),
                             ggarrange(ggplot() + thmHeatmapMeiS + 
                                         theme(legend.position='none') +
                                         ylab('scRT-Sim') + 
                                         xlab('') + 
                                         theme(axis.text.x=element_blank()),
                                       meiS$meanRT + 
                                         coord_cartesian(ylim=c(-2.5,2.5),xlim=c(0,129.4),expand=FALSE) +
                                         xlab('') + 
                                         ylab('normalized RT') + 
                                         scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c(-2,'',0,'',2)) +
                                         theme(plot.margin=unit(c(0,0,0.1,0.3),'cm'),
                                               axis.text.x=element_blank()),
                                       ncol=1,nrow=2,
                                       align='v',
                                       heights=c(1,1)))                           
gLeg <- getGGLegend(meiS$lstHM$fig + thmHeatmapMeiS)

## Heatmap: 
lstMetrics <- drawModelFitTable()

## J : S-phase split sort ---------------------------------------------------------
dAll <- lstMetrics$data

dAll$type                                   <- 'Others'
dAll$type[grep('MeiS',dAll$bg)]             <- 'Meiotic'
dAll$type[grep('Spg' ,dAll$bg)]             <- 'Gonia'
dAll$type[grep('(ESC|E14)' ,dAll$bg)]       <- 'ES Cells'
dAll$type[grep('(SSC|PGC)' ,dAll$bg)]       <- 'Germ Cells'
dAll$type[grep('(4C|2C|Sertoli)' ,dAll$bg)] <- 'Control'

dAll$type <- factor(dAll$type,
                    levels=c('Control','Meiotic','Gonia','ES Cells','Germ Cells','Others'))

lstData2 <- renameSamplesAndSplit(dAll)

gAll <- drawModelMetrics(lstData2$dfA[lstData2$dfA$topModels == 'Best models',],
                         myAlpha=.9)

gSort <- drawSortSections()

noMargin <- theme(plot.margin = margin(t=0, r=0, b=0, l=0, unit="cm"))
noXtxt   <- theme(axis.text.x=element_blank())

gT14data  <- ggarrange(gAll$gT14+ theme(plot.margin=margin(0,0,0,0,'cm')),
                       gAll$gT14S+ theme(plot.margin=margin(0,0,0,0,'cm'),
                                         axis.text.y=element_blank()),
                       ncol=2,nrow=1,
                       widths=c(4,4),
                       labels = c('',''),
                       align='h',
                       font.label = list(size=8),
                       hjust=0,vjust=1)

gSphaseSplit   <- ggarrange(gSort + theme(plot.margin=margin(0,0,0,0,'cm')),
                       gT14data + theme(plot.margin=margin(0,0,0,0,'cm')),
                       ncol=1,nrow=2,
                       widths=c(3,4),
                       labels = c('',''),
                       font.label = list(size=8),
                       hjust=0,vjust=1)

# ARRANGE PLOTS ################################################################
## A/C : RT + heatmaps ---------------------------------------------------------
gAC <- ggarrange(gLeg,meiS$gMeanSC,
                 esc$gMeanSC,
                 ncol=1,nrow=3,
                 heights=c(1,6,6),
                 labels=c('A','','C'),
                 font.label = list(size=8,face='bold'),
                 hjust=0,vjust=1)

## A/C : RT + NO heatmaps (allow PDF editing) ----------------------------------
gACn <- ggarrange(gLeg,meiS$gMeanSCemptyHM,
                 esc$gMeanSCemptyHM,
                 ncol=1,nrow=3,
                 heights=c(1,6,6),
                 labels=c('A','','C'),
                 font.label = list(size=8,face='bold'),
                 hjust=0,vjust=1)

## Panel B: Model fit heatmap --------------------------------------------------
gABC <- ggarrange(gAC,
                  lstMetrics$figs$gNoLbl,
                  ncol=2,nrow=1,
                  widths=c(16,6),
                  labels=c('','B'),
                  font.label = list(size=8,face='bold'),
                  hjust=0,vjust=1)

## JKL : Model props ridge plots------------------------------------------------
gDE   <- ggarrange(gAll$gRT+ theme(plot.margin=margin(0,0,0,0,'cm')),
                   gAll$gOri + theme(plot.margin=margin(0,0,0,0,'cm'),
                                     axis.text.y=element_blank()),
                   ncol=2,nrow=1,
                   widths=c(6,5),
                   labels = c('E','F'),
                   align='h',
                   font.label = list(size=8),
                   hjust=0,vjust=1)

gDEF   <- ggarrange(gSphaseSplit,
                    gDE,
                    ncol=2,nrow=1,
                    widths=c(1,2),
                    labels = c('D',''),
                    align='h',
                    font.label = list(size=8),
                    hjust=0,vjust=1)

gFinalPlot <- ggarrange(gABC,gDEF,
                        heights=c(1.3,1),
                        ncol=1,nrow=2)

# DRAW COMPLETE FIGURE 3 ######################################################
ggsave(filename = paste0('Pratto_et_al_Figure4.png'),
       plot = gFinalPlot, 
       height=6.5,
       width=6.5,
       dpi=400)

ggsave(filename = paste0('Pratto_et_al_Figure4.pdf'),
       plot = gFinalPlot, 
       height=6.5,
       width=6.5)

# DRAW COMPLETE FIGURE 3 (No heatmaps PDF) ######################################################
## Panel B: Model fit heatmap --------------------------------------------------
gABCn <- ggarrange(gACn,
                  lstMetrics$figs$gNoLbl,
                  ncol=2,nrow=1,
                  widths=c(16,6),
                  labels=c('','B'),
                  font.label = list(size=8,face='bold'),
                  hjust=0,vjust=1)

gFinalPlotN <- ggarrange(gABCn,gDEF,
                        heights=c(1.3,1),
                        ncol=1,nrow=2)
ggsave(filename = paste0('Pratto_et_al_Figure4_noHM.pdf'),
       plot = gFinalPlotN, 
       height=6.5,
       width=6.5)

# DRAW Supplementary Figure: Model properties ##################################

## Panel D-I: Model properties -------------------------------------------------
colz <- c('#555555','#de77ae')

# ---- Rep Rate
dfRR <- rbind(meiS$lstRepRate$data,esc$lstRepRate$data)

gSimTime <- ggplot(dfRR,aes(x=time,y=pc,color=type)) + 
  geom_line(lwd=.3) + 
  geom_point(size=.3) + 
  geom_hline(yintercept=90,lwd=.2,color='magenta',lty='dashed') + 
  ylab('Replicated (%)') + 
  xlab('Time (Simulation cycles)') + 
  annotate(geom='text',
           label=paste0('90%'),
           size=7*5/14,
           x=10,y=90,
           hjust=0,vjust=1.4,
           color='magenta',
           check_overlap=TRUE) + 
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) +
  theme(legend.position=c(1,0),
        legend.justification=c(1,0),
        legend.background=element_blank(),
        legend.key.size=unit(0.3,'cm'),
        legend.title=element_blank(),
        legend.text = element_text(size=7))

# ---- Forks used per cycle
dfForks <- rbind(meiS$lstForksUsed$data, esc$lstForksUsed$data)

gActiveForks <- ggplot(dfForks,aes(x=as.numeric(time),y=nf,color=type)) + 
  geom_point(alpha=1,size=.2) +
  ylab('Active forks (#)') + 
  xlab('Time (Model cycles)') + 
  theme(legend.position='none')+
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) 

# ---- Origins used per S-phase
dfOriFiringFreq <- rbind(meiS$lstOriPerSim$data, esc$lstOriPerSim$data)

gTotOri <- ggplot(dfOriFiringFreq,aes(x=nOris,fill=type)) + 
  geom_histogram(binwidth = 25,color='grey40',lwd=.3) + 
  xlab('Origins fired (#)') + 
  ylab('Simulations (%)') + 
  geom_label(data = dfOriFiringFreq %>% group_by(type) %>% 
               summarize(med=round(median(nOris))) %>% 
               mutate(row = row_number()),
             aes(x=mean(dfOriFiringFreq$nOris)-(mean(dfOriFiringFreq$nOris)-med)/3,
                 y=2+(row-1)*10,
                 label=paste0('Median:\n',format(med,big.mark = ",")," origins"),
                 color=type),
             fill='white',
             size=7*5/14,
             hjust=.5,vjust=0,
             label.size=NA) +
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) +
  theme(legend.position='none')

# ---- Distribution of origin firing freq.s
dfOriFreqDist <- rbind(meiS$lstOriFiring$data, esc$lstOriFiring$data)

gOriFreq <- ggplot(dfOriFreqDist) + 
  geom_line(aes(x=freq,y=pc,group=type,color=type),lwd=.3) + 
  geom_point(aes(x=freq,y=pc,group=type,color=type),size=.3) + 
  ylab('Origins (%)') + 
  xlab('Firing frequency\n(Percent of simulations)') + 
  geom_vline(data=dfOriFreqDist %>% group_by(type) %>% 
               summarize(med=mean(medianFreq)) %>% mutate(row=row_number()) ,
             aes(xintercept=med,color=type),
             lwd=.2,lty='dashed') + 
  geom_label(data=dfOriFreqDist %>% group_by(type) %>% 
               summarize(med=mean(medianFreq)) %>% mutate(row=row_number()) ,
             aes(x=med,
                 label=paste0('Median = ',med,' %'),
                 color=type,
                 y=max(dfOriFreqDist$pc)*(10-(row*2)+2)/10),
             size=7*5/14,
             hjust=-0.1,vjust=1,
             label.size=NA) + 
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) +
  theme(legend.position='none')

# ---- Fork travel distance
dfForkDist <- rbind(meiS$lstOriUsage$data, esc$lstOriUsage$data)

gForkDist <- ggplot(dfForkDist) + 
  geom_line(aes(x=forkDist,y=pc,color=type),lwd=.3) + 
  geom_point(aes(x=forkDist,y=pc,color=type),size=.3) + 
  geom_vline(data=dfForkDist %>% group_by(type) %>% 
               summarize(med=mean(medianDist)) %>% mutate(row=row_number()) ,
             aes(xintercept=med,color=type),
             lwd=.2,lty='dashed') + 
  geom_label(data=dfForkDist %>% group_by(type) %>% 
               summarize(med=mean(medianDist)) %>% mutate(row=row_number()) ,
             aes(x=med,
                 label=paste0('Median = ',med,' Mb'),
                 color=type,
                 y=max(dfForkDist$pc)*(10-(row*2)+2)/10),
             size=7*5/14,
             hjust=-0.1,vjust=1,
             label.size=NA) + 
  ylab('Forks (%)') + 
  xlab('Distance replicated (Mb)') + 
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) +
  theme(legend.position='none')

# ---- New ori firing rate per cycle
dfNewOris <- rbind(meiS$lstNewOriRate$data, esc$lstNewOriRate$data)

gOriPerCycle <- ggplot(dfNewOris) + 
  scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),
                labels=c(0.01,0.1,1,10,100,1000)) + 
  geom_point(aes(x=vTimeFired,y = n/100,color=type,group=type),size=.3) + 
  geom_line(aes(x=vTimeFired,y = n/100,color=type,group=type),lwd=.3) + 
  geom_hline(data=dfNewOris %>% group_by(type) %>% 
               summarize(med=median(n)/100) %>% mutate(row=row_number()) ,
             aes(yintercept=med,
                 color=type),
             lty='dashed',
             lwd=.3) + 
  xlab('Time fired (cycle)') + 
  annotation_logticks(sides='l',
                      long  = unit(0.1,'cm'),
                      mid   = unit(0.05,'cm'),
                      short = unit(0.05,'cm'),
                      size=.2) +
  ylab('New origins fired\n(Mean count per simulation)') +
  geom_label(data=dfNewOris %>% group_by(type) %>% 
               summarize(med=median(n)/100) %>% mutate(row=row_number()),
             x=Inf,
             aes(y=10^(row+1),
                 color=type,
                 label=paste0('Median = ',round(med,1),' origins')),
             size=7*5/14,
             hjust=1,vjust=1,
             label.size=NA)+
  scale_fill_manual(values=colz) + 
  scale_color_manual(values=colz) +
  theme(legend.position='none')

# ---- All together now !
gSAll <- ggarrange(gSimTime,gTotOri,gOriFreq,gActiveForks,gForkDist,gOriPerCycle,
                 ncol=3,nrow=2,
                 align='v',
                 labels=c('B','D','F','C','E','G'),
                 font.label = list(size=8,face='bold'),
                 hjust=0,vjust=1)

gRTEML <- drawRTEMLFig()

gSuppProps <- ggarrange(gRTEML$fig, gSAll,
                        ncol=1,nrow=2,
                        labels=c('A',''),
                        font.label = list(size=8,face='bold'),
                        hjust=0,vjust=1)

ggsave('Pratto_et_al_SupplementalModelProperties_forFigure4.png',gSuppProps,width=6.5,height=8,dpi = 300)
ggsave('Pratto_et_al_SupplementalModelProperties_forFigure4.pdf',gSuppProps,width=6.5,height=8)

# DRAW Supplementary Figure: Ridgeplots with SPO11-ko ##########################
gSuppX2 <- ggarrange(esc$genomeRT + ggtitle('ESC S-phase (a) : Experimental vs simulated RT'),
                     meiS$genomeRT + ggtitle('Meiotic S-phase (a) : Experimental vs simulated RT'),
                     ncol=1,nrow=2,
                     heights=c(3,3),
                     labels = c('A','B'),
                     font.label = list(size=8),
                     hjust=0,vjust=1)

ggsave('Pratto_et_al_Supplemental_RTsim_vs_RTexp.png',gSuppX2,height=9,width=6.5,dpi=300)
ggsave('Pratto_et_al_Supplemental_RTsim_vs_RTexp.pdf',gSuppX2,height=9,width=6.5)

# DRAW Supplementary Figure: Ridgeplots with SPO11-ko ##########################
gTop  <- drawSpo11KOSortingSchema(.1,.1)

gLeft <- ggarrange(gAll$gTable + theme(plot.margin=margin(0,0,0,0,'cm')),
                   gAll$gRTrel + theme(plot.margin=margin(0,0,0,0,'cm')),
                   ncol=1,nrow=2,
                   heights=c(3,3),
                   labels = c('C','D'),
                   font.label = list(size=8),
                   hjust=0,vjust=1)

gRight <- ggarrange(gAll$gSR2  + theme(plot.margin=margin(0,0,0,0,'cm')),
                    gAll$gSRT  + theme(plot.margin=margin(0,0,0,0,'cm')),
                    gAll$gSOri + theme(plot.margin=margin(0,0,0,0,'cm')),
                    ncol=1,nrow=3,
                    labels = c('E','F','G'),
                    font.label = list(size=8),
                    hjust=0,vjust=1)

gSupplementLR <- ggarrange(gLeft,gRight,widths=c(2.2,2),ncol=2,nrow=1)
gSupplement   <- ggarrange(gTop,gSupplementLR,heights=c(1.5,4),ncol=1,nrow=2)

ggsave('Pratto_et_al_Spo11SupplementalRidgePlot_forFigure4.png',width=6.5,height=9,dpi = 300)
ggsave('Pratto_et_al_Spo11SupplementalRidgePlot_forFigure4.pdf',width=6.5,height=9)

###
write.csv(gAll$ssMain,file  = 'summarystatsForMainRidgePlot.csv',col.names = TRUE,row.names = FALSE)
write.csv(gAll$ssSpo11,file = 'summarystatsForSpo11RidgePlot.csv',col.names = TRUE,row.names = FALSE)


