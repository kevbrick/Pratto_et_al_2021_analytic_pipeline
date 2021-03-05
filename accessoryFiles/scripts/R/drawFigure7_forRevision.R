imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')
source('accessoryFiles/scripts/R/repliSim_loadModules.R')

library(data.table)

theme7point()

# FUNCTIONS --------------------------------------------------------------------
getHumanSort1 <- function(fcsIn ='Specimen_001_Human-3_004.fcs', limz = NULL) {
  dSort        <- importFCS(fcsIn)
  
  names(dSort)[names(dSort) == 'BV421.A'] <- 'DAPI'
  names(dSort)[names(dSort) == 'B530.A'] <- 'HORMAD1'
  names(dSort)[names(dSort) == 'R660.A'] <- 'DMRT1'
  
  
  dSort$DAPI2C  <- dSort$DAPI  > 90000 & dSort$DAPI < 125000
  dSort$DAPIRep  <- dSort$DAPI  > 110000 & dSort$DAPI < 140000
  #dSort$DAPIstrRep  <- dSort$DAPI  > 110000 & dSort$DAPI < 120000
  dSort$DAPI4C  <- dSort$DAPI  > 150000 & dSort$DAPI < 220000
  dSort$DAPIok  <- dSort$DAPI  > 30000 & dSort$DAPI < 220000 & dSort$BV421.W < 100000
  dSort$H1ok <- dSort$HORMAD1 > 1500
  dSort$scatterOK <- dSort$FSC.A < 2e5 & dSort$SSC.A < 2e5
  
  dSort$DSok    <- dSort$DAPIok 
  
  dSort$allOK             <- 'Other'
  dSort$allOK[dSort$DSok] <- 'DAPI_Gated'
  dSort$ploidy <- '1C'
  dSort$ploidy[dSort$DAPI2C] <- '2C'
  dSort$ploidy[dSort$DAPIRep] <- 'Rep'
  #dSort$ploidy[dSort$DAPIstrRep] <- 'StrRep'
  dSort$ploidy[dSort$DAPI4C] <- '4C'
  
  #dSortOK <- dSort[dSort$DSok,]
  dSortOK <- dSort[dSort$DSok & dSort$scatterOK,]
  return(list(all = dSort, ok=dSortOK))
  
}

getHumanSort2 <- function(fcsIn ='2019-12-04_Human_002.fcs', limz = NULL) {
  
  dSort        <- importFCS(fcsIn)
  names(dSort)[names(dSort) == 'BV421.A'] <- 'DAPI'
  names(dSort)[names(dSort) == 'B530.A'] <- 'SYCP3'
  names(dSort)[names(dSort) == 'R660.A'] <- 'DMRT1'
  
  
  dSort$DAPI2C  <- dSort$DAPI  > 90000 & dSort$DAPI < 120000
  dSort$DAPIRep  <- dSort$DAPI  > 120000 & dSort$DAPI < 160000
  dSort$DAPI4C  <- dSort$DAPI  > 160000 & dSort$DAPI < 220000
  dSort$DAPIok  <- dSort$DAPI  > 30000 & dSort$DAPI < 220000 & dSort$BV421.W < 80000
  
  dSort$scatterOK <- dSort$FSC.A < 2e5 & dSort$SSC.A < 2e5
  dSort$DSok    <- dSort$DAPIok 
  
  dSort$allOK             <- 'Other'
  dSort$allOK[dSort$DSok] <- 'DAPI_Gated'
  dSort$ploidy <- '1C'
  dSort$ploidy[dSort$DAPI2C] <- '2C'
  dSort$ploidy[dSort$DAPIRep] <- 'Rep'
  dSort$ploidy[dSort$DAPI4C] <- '4C'
  
  dSortOK <- dSort[dSort$DSok & dSort$scatterOK,]
  
  return(list(all = dSort, ok=dSortOK))
  
}

plotSort1     <- function(colz    = c('salmon','dodgerblue3','salmon','grey80'), pointSz = 0.2, pointA  = 0.5){
  
  theme7point()
  
  sd1 <- getHumanSort1()
  
  dRep <- sd1$ok[sd1$ok$DAPIRep,]
  d4C <- sd1$ok[sd1$ok$DAPI4C,]
  d2C <- sd1$ok[sd1$ok$DAPI2C,]
  nTot <- length(sd1$ok$DAPI)
  nW <- 750
  
  dDAPI    <- as.data.frame(table(round(sd1$ok$DAPI/nW)))
  dDAPIRep <- as.data.frame(table(round(sd1$ok$DAPI[sd1$ok$DAPIRep]/nW)))
  
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
  
  dLabels <- data.frame(x=c(5,9.3,13,17.5),
                        y=c(0.05,0.05,1.0,0.05),
                        lbl=c('1C','2C','2-4C','4C'))
  
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
  
  
  dRep$stage <- "NA"
  dRep$stage[dRep$DMRT1 > 8000 & dRep$DMRT1 < 60000 & dRep$HORMAD1 < 500  & dRep$HORMAD1 > 20 ] <- 'Spg'
  dRep$stage[dRep$DMRT1 > 1000 & dRep$DMRT1 < 8000 & dRep$HORMAD1 < 500  & dRep$HORMAD1 > 20] <- 'MeiS'
  dRep$stage[dRep$DMRT1 < 300 & dRep$DMRT1 > 10 & dRep$HORMAD1 < 2000  & dRep$HORMAD1 > 300] <- 'LMeiS_PreL'
  
  dRep$stage <- factor(dRep$stage, levels=c('Spg','MeiS','LMeiS_PreL','NA'))
  
  gSort <- ggplot(dRep,aes(x=HORMAD1,y=DMRT1,
                           fill=as.factor(stage),
                           color=as.factor(stage))) + 
    geom_density2d(lwd=.05,alpha=.8) +
    geom_point(shape=20,size=pointSz,alpha=pointA) +
    scale_x_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) + 
    scale_y_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) +  
    annotate(geom='rect',
             xmax=500,xmin=20,
             ymax=60000,ymin=9000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[1])+
    annotate(geom='rect',
             xmax=500,xmin=20,
             ymax=8000,ymin=1000,
             lwd=.3, lty='solid',
             fill=alpha('white',0),color=colz[2])+
    annotate(geom='rect',
             xmax=2000,xmin=300,
             ymax=300,ymin=10,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[3])+
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
    xlab('HORMAD1 fluorescence') +
    ylab('DMRT1 fluorescence') +
    coord_cartesian(ylim=c(10,110000),xlim=c(10,3000))+
    annotate(geom='text',
             x=20,y=100000,
             label='Spg',
             hjust=0,
             size=7*5/14,
             color=colz[1]) +
    annotate(geom='text',
             x=20,y=800,
             label='MeiS & Spg',
             hjust=0,vjust=1.5,
             size=7*5/14,
             color=colz[2])+ 
    annotate(geom='text',x=270,y=15,
             label='MeiP / Late MeiS',
             hjust=1,vjust=0,
             size=7*5/14,
             color=colz[3]) 
  
  return(list(gMain = gSort,
              gDAPI = gDNAContent))
}

plotSort2 <- function(colz    = c('darkolivegreen4','dodgerblue3','grey'), pointSz = 0.2, pointA  = 0.5){
  
  theme7point()
  
  sd1 <- getHumanSort2()
  
  dRep <- sd1$ok[sd1$ok$DAPIRep,]
  d4C <- sd1$ok[sd1$ok$DAPI4C,]
  d2C <- sd1$ok[sd1$ok$DAPI2C,]
  nTot <- length(sd1$ok$DAPI)
  nW <- 750
  
  dDAPI    <- as.data.frame(table(round(sd1$ok$DAPI/nW)))
  dDAPIRep <- as.data.frame(table(round(sd1$ok$DAPI[sd1$ok$DAPIRep]/nW)))
  
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
  
  dLabels <- data.frame(x=c(5,9.3,13,17.5),
                        y=c(0.05,0.05,1.0,0.05),
                        lbl=c('1C','2C','2-4C','4C'))
  
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
  
  
  dRep$stage <- "NA"
  #dRep$stage[dRep$DMRT1 > 8000 & dRep$DMRT1 < 60000 & dRep$SYCP3 < 200  & dRep$SYCP3 >20  ] <- 'p1'
  #dRep$stage[dRep$DMRT1 > 2000 & dRep$DMRT1 < 8000 & dRep$SYCP3 < 500  & dRep$SYCP3 > 200] <- 'p2'
  dRep$stage[dRep$DMRT1 > 200 & dRep$DMRT1 < 2000 & dRep$SYCP3 < 500  & dRep$SYCP3 > 200] <- 'MeiS & Spg'
  dRep$stage[dRep$DMRT1 < 200 & dRep$DMRT1 > 0 & dRep$SYCP3 < 500  & dRep$SYCP3 > 200] <- 'MeiS'
  
  dRep$stage <- factor(dRep$stage, levels=c('MeiS','MeiS & Spg',"NA"))
  
  gSort <- ggplot(dRep,aes(x=SYCP3,y=DMRT1,
                           fill=as.factor(stage),
                           color=as.factor(stage))) + 
    geom_density2d(lwd=.05,alpha=.8) +
    geom_point(shape=20,size=pointSz,alpha=pointA) +
    scale_x_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) + 
    scale_y_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) +  
    annotate(geom='rect',
             xmax=500,xmin=200,
             ymax=2000,ymin=200,
             lwd=.3, lty='solid',
             fill=alpha('white',0),color=colz[2])+
    annotate(geom='rect',
             xmax=500,xmin=200,
             ymax=200,ymin=20,
             lwd=.3, lty='solid',
             fill=alpha('white',0),color=colz[1])+
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
    xlab('SYCP3 fluorescence') +
    ylab('DMRT1 fluorescence') +
    coord_cartesian(ylim=c(10,110000),xlim=c(30,1000))+
    annotate(geom='text',
             x=180,y=400,
             label='MeiS & Spg',
             hjust=1,vjust=0,
             size=7*5/14,
             color=colz[2]) +
    annotate(geom='text',
             x=180,y=40,
             label='MeiS',
             hjust=1,vjust=0,
             size=7*5/14,
             color=colz[1])
  
  return(list(gMain = gSort,
              gDAPI = gDNAContent))
}

plotSCExpression <- function(){
  dfSCRNA          <- read.table('scRNASeq_Human.tab',header=TRUE)
  
  dfSCRNA$cellType <- factor(dfSCRNA$cellType,
                             levels=c('Spg(Undiff)','Spg(Diff)','Meiotic(Early)','Meiotic(Late)'))
  
  gExp <- ggplot(dfSCRNA,
                 aes(x=cellType,y=meanExpression*percentCells/100,color=gene,group=gene)) + 
    geom_line(lwd=.4) + 
    geom_point(aes(shape=gene),fill='white') + 
    scale_shape_manual(values=c(21,22,24,25)) +
    xlab('') + 
    ylab('scRNA-Seq\nexpression level') + 
    theme(legend.position=c(1,.5),
          legend.title=element_blank(),
          legend.justification=c(1,.5),
          legend.key.size=unit(0.2,'cm'),
          axis.text.x=element_text(angle=25,hjust=1))
  
  return(gExp)
}

drawD2T <- function(df,lst,isMain=TRUE){
  dfMlt <- reshape2:::melt.data.frame(df,id.vars='d2tel', 
                                      measure.vars = lst,
                                      variable.name = "celltype",
                                      value.name = "RTn") %>%
    mutate(d2rnd = floor(d2tel/2)*2) 
  
  dfPlot <- as.data.frame(dfMlt %>% group_by(celltype,d2rnd) %>% 
                            summarise(mean = mean(RTn), 
                                      median=median(RTn), 
                                      ln = mean(RTn) - sd(RTn), 
                                      hn = mean(RTn) + sd(RTn)))
  
  gInit <- ggplot(dfPlot %>% dplyr::filter(d2rnd < 51),
                 aes(x=d2rnd,
                     fill=celltype,
                     color=celltype)) + 
    geom_ribbon(aes(ymin=ln,ymax=hn,alpha=.1),lwd=.2) +  
    geom_boxplot(data=dfMlt %>% dplyr::filter(d2rnd < 51),
                 aes(x=d2rnd,y=RTn,group=d2rnd),
                 color='grey30',
                 lwd=.2,
                 outlier.size=.01,
                 outlier.alpha=0) + 
    geom_point(aes(y=mean),size=.1,color='firebrick') + 
    theme(legend.position='none',
          strip.text=element_text(hjust=0,size=7),
          strip.background=element_blank()) + 
    coord_cartesian(ylim=c(0,1),xlim=c(-2,52), expand=FALSE) + 
    geom_hline(yintercept=0,lwd=.2,color='black') + 
    xlab('Distance to telomere (Mb)') + 
    ylab('Normalized RT')
  
  gInit <- ggplot(dfPlot %>% dplyr::filter(d2rnd < 51),
                  aes(x=d2rnd,
                      fill=celltype,
                      color=celltype)) + 
    geom_ribbon(aes(ymin=ln,ymax=hn,alpha=.1),lwd=.2) +  
    geom_boxplot(data=dfMlt %>% dplyr::filter(d2rnd < 51),
                 aes(x=d2rnd,y=RTn,group=d2rnd),
                 color='grey30',
                 lwd=.2,
                 outlier.size=.01,
                 outlier.alpha=0) + 
    geom_point(aes(y=mean),size=.3,color='magenta') + 
    theme(legend.position='none',
          strip.text=element_text(hjust=0,size=7),
          strip.background=element_blank()) + 
    coord_cartesian(ylim=c(0,1),xlim=c(-2,52), expand=FALSE) + 
    geom_hline(yintercept=0,lwd=.2,color='black') + 
    xlab('Distance to telomere (Mb)') + 
    ylab('Normalized RT')

  if (isMain){
    gRet <- gInit + facet_wrap(~celltype,ncol=2) + 
      scale_fill_manual(values=c('darkolivegreen3',rep('dodgerblue2',3),rep('grey80',6))) + 
      scale_color_manual(values=c('darkolivegreen4',rep('dodgerblue3',3),rep('grey50',6))) 
  }else{
    gRet <- gInit + facet_wrap(~celltype,ncol=2) + 
      scale_fill_manual(values=c('darkolivegreen3',rep('dodgerblue2',3),rep('grey80',6))) + 
      scale_color_manual(values=c('darkolivegreen4',rep('dodgerblue3',3),rep('grey50',6))) 
  }  
  return(gRet)
}

drawGC  <- function(df){
  dfMlt <- reshape2:::melt.data.frame(df,id.vars='d2tel', measure.vars = 'GC-content')
  
  dfMlt$d2rnd <- floor(dfMlt$d2tel/2)*2
  
  dfPlot <- as.data.frame(dfMlt %>% group_by(variable,d2rnd) %>% 
                            summarise(mean = mean(value), 
                                      median=median(value), 
                                      ln = mean(value) - sd(value), 
                                      hn = mean(value) + sd(value)))
  
  gRet <- ggplot(dfPlot[dfPlot$d2rnd < 51,],
                  aes(x=d2rnd,fill=variable,color=variable)) + 
    geom_ribbon(aes(ymin=ln,ymax=hn,alpha=.1),lwd=.2) +  
    geom_boxplot(data=dfMlt[dfMlt$d2rnd < 51,],
                 aes(x=d2rnd,y=value,group=d2rnd),
                 color='grey30',
                 lwd=.2,
                 outlier.size=.01,
                 outlier.alpha=0) + 
    geom_point(aes(y=mean),size=.1,color='firebrick') + 
    theme(legend.position='none',
          strip.text=element_text(hjust=0,size=7),
          strip.background=element_blank()) + 
    coord_cartesian(ylim=c(25,75),xlim=c(-2,52), expand=FALSE) + 
    geom_hline(yintercept=0,lwd=.2,color='black') + 
    xlab('Distance to telomere (Mb)') + 
    ylab('GC-content (%)') + 
    scale_fill_manual(values=c('pink')) + 
    scale_color_manual(values=c('grey40')) 
  
  return(gRet)
}

normEm <- function(x,fld,nm){
  
  xMin <- min(x[[fld]])
  xMax <- max(x[[fld]])
  
  x$oneTo100 <- (x[[fld]] - min(x[[fld]])) / max((x[[fld]] - min(x[[fld]])))
  print(min(x$oneTo100))
  print(max(x$oneTo100))
  x        <- x[x$d2tel <= 45,]
  x$std    <- standardizeMNSD(x[[fld]])
  x$std    <- x[[fld]]
  x$std    <- x$oneTo100
  
  #x        <- x[x[[fld]] > quantile(x[[fld]],0.05) & x[[fld]] < quantile(x[[fld]],0.95)]
  gc       <- aggregate(x$pcGC  ,by=list(d = x$d2tel),FUN=mean)
  i        <- aggregate(x$std   ,by=list(d = x$d2tel),FUN=mean)
  iSD      <- aggregate(x$std   ,by=list(d = x$d2tel),FUN=sd)
  iMean   <- aggregate(x$std   ,by=list(d = x$d2tel),FUN=mean)
  iSum     <- aggregate(x$std   ,by=list(d = x$d2tel),FUN=sum)
  iMed     <- aggregate(x$std   ,by=list(d = x$d2tel),FUN=median)
  
  i$med    <- iMed$x
  i$mean   <- iMean$x
  i$sd     <- iSD$x
  i$sum    <- iSum$x
  i$nm     <- nm

  i$mnsd   <- standardizeMNSD(i$x)
  i$mnsdsd <- standardizeMNSD(i$sd)
  i$mnsdGC <- standardizeMNSD(i$x/gc$x)

  i$fc     <- abs(i$x)/mean(abs(i$x))
  if (grepl('simRT',fld)){
    i$fc <- 1/i$fc
  }

  return(i)
}

drawCCplot <- function(df){

  ijAll         <- rbind(drawCC_v_RR(df,'SSDSaa1','SSDScl4','A/A v C/L4'),
                         drawCC_v_RR(df,'SSDSaa1','SSDSaa2','A/A v A/A'))

  ijAll$species <- gsub('_.+'     ,'', toupper(ijAll$name))
  ijAll$species <- gsub('MM'     ,'Mouse', ijAll$species)
  ijAll$species <- gsub('HS'     ,'Human', ijAll$species)
  ijAll$comp    <- gsub('(MM|HS)_','', toupper(ijAll$name))
  ijAll$comp    <- gsub('V',' V ', toupper(ijAll$comp))

  gCC <- ggplot(ijAll,aes(x=win*10000/1000000,y=cc12,
                          color=name,fill=name,linetype=name)) +
    geom_point(size=.6) + geom_line(lwd=.3) +
    geom_vline(xintercept=1,lty='dashed',lwd=.2) +
    scale_x_log10() +
    scale_y_continuous(breaks=seq(0,1,by=.2),
                       labels=c(0,'','','','',1),
                       position='left') +
    annotation_logticks(sides='b',
                        long=unit(0.2,'cm'),
                        mid=unit(0.1,'cm'),
                        short=unit(0.1,'cm'),
                        size=.2) +
    xlab('Window size (Mb)') +
    ylab(bquote("Correlation between DSB maps"*" ("*R^2*")")) +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background=element_blank(),
          legend.title=element_blank(),
          strip.background=element_blank(),
          strip.text=element_text(size=8,face='bold')) +
    scale_color_manual(values=c('red','red')) +
    scale_fill_manual(values=c('red','red')) +
    scale_linetype_manual(values=c(1,2)) +
    coord_cartesian(ylim=c(0,1),xlim=c(0.008,22),expand=FALSE)

  return(gCC)
}

drawGCCC   <- function(df){

  hColz      <- names(df)[grep('exp',names(df))]
  hS         <- apply(df[,..hColz],1,function(x){sum(abs(x) < 0.001)>1})
  hDF        <- as.data.frame(df[!hS,])

  ccH                                          <- as.data.frame(cor(hDF[,c(grep('expRT|pcGC',names(hDF)))]))
  ccHS                                         <- data.frame(nm=rownames(ccH),R=ccH$pcGC)
  ccHS                                         <- ccHS[!(ccHS$nm %in% c("pcGC","expRT_MeiS_weakH1_hsS1_2to4C_HORMAD1_DMRT1_hg38")),]
  ccHS$type                                    <- 'Others'
  ccHS$type[grep('SPG|Spg|Mei|Germ',ccHS$nm)]  <- 'Germ-line'
  ccHS$type[grep('_(2C|4C|Sertoli)_',ccHS$nm)] <- 'Non-Rep.'
  ccHS$species                                 <- 'Human'

  ccHS <- ccHS[ccHS$type != 'Non-Rep.',]
  ccAll      <- rbind(ccHS)

  ccAll$type <- factor(ccAll$type,levels=c('Germ-line','Others','Non-Rep.'))

  ### STATS

  oMH <- ccAll$R[ccAll$type=='Others' & ccAll$species == 'Human']
  gMH <- ccAll$R[ccAll$type=='Germ-line' & ccAll$species == 'Human']

  wH <- wilcox.test(gMH,oMH)

  toSci <- function(x){
    xChr <- as.character(x)
    if(grepl("e",xChr)){
      gsub()
    }else{
      return(x)
    }
  }

  ccAll$Pval <- "P = ??"
  #ccAll$Pval[ccAll$species == 'Mouse'] <- paste0('P = ',sprintf("%.4f",wM$p.value))
  ccAll$Pval[ccAll$species == 'Human'] <- paste0('P = ',sprintf("%.4f",wH$p.value))

  gGCCC <- ggplot(ccAll,aes(x=type,y=R,fill=type)) +
    geom_boxplot(lwd=.2,outlier.size=0,outlier.alpha=0,varwidth=TRUE,alpha=.4) +
    geom_point(size=1,shape=21,lwd=.1,aes(color=type),
               fill='black',position=position_jitter(seed = 12, height = 0,width=.1)) +
    theme(legend.position='none',
          strip.background=element_blank(),
          strip.text=element_text(size=8,face = 'bold')) +
    xlab('') + ylab('Correlation\nwith GC content\n(Pearson R)') +
    geom_text(aes(label=Pval),
              x=1.5,y=0.8,check_overlap = TRUE,
              hjust=0.5,size=7*5/14) +
    scale_fill_manual(values=c('forestgreen','grey50')) +
    scale_color_manual(values=c('forestgreen','grey50')) +
    annotate(geom='segment',x=1  ,xend=2  ,y=0.75, yend=0.75  ,lwd=.2) +
    annotate(geom='segment',x=1  ,xend=1  ,y=0.75, yend=0.7   ,lwd=.2) +
    annotate(geom='segment',x=1.5,xend=1.5,y=0.75, yend=0.78  ,lwd=.2) +
    annotate(geom='segment',x=2  ,xend=2  ,y=0.75, yend=0.7   ,lwd=.2)  +
    coord_cartesian(ylim=c(0.35,0.85))
    scale_y_continuous(breaks=c(0.4,0.6,0.8))

  return(gGCCC)
}

# Make Sorting Figs ------------------------------------------------------------
gPanelSCexpression       <-  plotSCExpression()

lstHumanSort1 <- plotSort1()
gPanelSort1       <- lstHumanSort1$gMain

lstHumanSort2 <- plotSort2()
gPanelSort2       <- lstHumanSort2$gMain

# Make other figs --------------------------------------------------------------
tInit      <- fread('rep_v_rec_HG38.tab',header=TRUE)
tColz      <- names(tInit)[grep('exp',names(tInit))]
tS         <- apply(tInit[,..tColz],1,function(x){sum(abs(x) < 0.001)>1})
t          <- tInit[!tS,]

### Get d2tel 
nSz <- 1000000
t$minD <- apply(t[,c('pDist','qDist')],1,min)
t$d2tel <- floor(t$minD/nSz)/(1000000/nSz)

# Draw CC plot -----------------------------------------------------------------
gPanelCC     <- drawCCplot(t)

# Draw DSBs vs Tels Plot -------------------------------------------------------
tAA1  <- normEm(t,'SSDSaa1','AA1')
tAA2  <- normEm(t,'SSDSaa2','AA2')
tCL4  <- normEm(t,'SSDScl4','CL4')

dfDSBs      <- rbind(tAA1,tAA2,tCL4)
dfDSBs$nm[dfDSBs$nm == 'AA1']     <- 'DSBs (A/A; 1)'
dfDSBs$nm[dfDSBs$nm == 'AA2']     <- 'DSBs (A/A; 2)'
dfDSBs$nm[dfDSBs$nm == 'CL4']     <- 'DSBs (C/L4)'

dfDSBs$nm <- factor(dfDSBs$nm,levels=c('DSBs (A/A; 1)','DSBs (A/A; 2)','DSBs (C/L4)'))

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colorBlindBlack3  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gPanelDSBTel <- ggplot(dfDSBs[dfDSBs$d <= 50,],
                  aes(x=d,y=fc,
                      color=nm,
                      fill=nm),
                  shape=21) +
  geom_vline(xintercept=10,color='grey50',lwd=.3,lty='dashed') +
  geom_point(size=1,lwd=.3) +
  geom_smooth(alpha=.1,span=.4,lwd=.4,se=FALSE) +
  ylab('Fold enrichment') +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.key.size=unit(0.3,'cm'),
        legend.title = element_blank()) +
  scale_linetype_manual(values=c(1,1,2,1,1)) +
  scale_color_manual(values=c('grey40','red','forestgreen','black','grey70')) +
  scale_fill_manual(values=c('grey20','firebrick','darkolivegreen3','black','grey70')) +
  scale_color_manual(values=colorBlindBlack3) +
  scale_fill_manual(values=colorBlindBlack3) +
  xlab('Distance to telomere (Mb)') +
  coord_cartesian(xlim=c(-1,46),expand=FALSE) +
  scale_y_continuous(position = 'left')

# Make RT vs position fig ------------------------------------------------------ 
t$'Meiotic'                         <- standardize0to1(2^t$expRT_MeiS_ALL_hsS1_2to4C_SCP3_NA)
t$'Meiotic & Spermatogonia (a)'     <- standardize0to1(2^t$expRT_Germline_ALL_S2_2to4C_SCP3_DMRT1w)
t$'Meiotic & Spermatogonia (b)'     <- standardize0to1(2^t$expRT_MeiS_ALL_hsS1_2to4C_HORMAD1DMRT1_NA)
t$'Meiotic & Spermatogonia (c)'     <- standardize0to1(2^t$expRT_MeiS_ALL_S2_2to4C_SCP3_DMRT1)
t$'RPE1 cells'                      <- standardize0to1(2^t$expRT_RPE1_MID1_RT)
t$'ES cells'                        <- standardize0to1(t$expRT_H7_ESC_Ext35479608)
t$'Myoblast cells'                  <- standardize0to1(t$'expRT_FM01-154-001_Myoblast_Int58331187')
t$'LCL cells'                       <- standardize0to1(t$expRT_GM06990_Lymphoblastoid_Ext54054609)
t$'Liver cells'                     <- standardize0to1(t$expRT_CyT49_Liver_D16_Int81158282)
t$'U2OS cells'                      <- standardize0to1(t$expRT_U2OS_Bone_Epithelial_Int66343918)
t$'GC-content'                      <- t$pcGC*100

gVSTel <- drawD2T(t,c('Meiotic',
                      'Meiotic & Spermatogonia (a)',
                      'Meiotic & Spermatogonia (b)',
                      'Meiotic & Spermatogonia (c)',
                      'RPE1 cells',
                      'ES cells',
                      'Liver cells',
                      'Myoblast cells',
                      'LCL cells',
                      'U2OS cells'),isMain = TRUE)

# ## Panel C:
# gMainD2T <- drawD2T(t,c('Meiotic','RPE1 cells'))
# gGC <- drawGC(t)
#
# gPanelC  <- ggarrange(gMainD2T,gGC,align='v',ncol=1,nrow=2,
#                       heights=c(2,1),
#                       labels=c('C',''),
#                       font.label = list(size=8),
#                       hjust=0,vjust=1)

# Make GCv Tel Fig ------------------------------------------------------------- 
gPanelG <- drawGC(t)

# Make GC boxplot Fig ---------------------------------------------------------- 
gPanelH <- drawGCCC(tInit) 

# Merge into figure ------------------------------------------------------------ 
ltMarg <- theme(plot.margin=margin(.2,0,0,.2,'cm'))
lMarg <- theme(plot.margin=margin(0,0,0,0.3,'cm'))
noMarg <- theme(plot.margin=margin(0,0,0,0,'cm'))

gPanelAB <- ggarrange(gPanelCC + lMarg,
                      gPanelDSBTel + noMarg,
                      labels = c("A","B"), 
                      ncol=2,nrow=1,
                      font.label = list(size=8),
                      hjust=0,vjust=1)

gPanelSort <- ggarrange(gPanelSCexpression + lMarg,
                       gPanelSort1 + noMarg,
                       gPanelSort2 + noMarg,
                       gPanelG + ltMarg,
                       gPanelH + ltMarg,
                       labels = c("C","D","E","G","H"), 
                       ncol=1,nrow=5,
                       heights=c(.8,1,1,.9,.9),
                       font.label = list(size=8),
                       hjust=0,vjust=1)


gPanelABF <- ggarrange(gPanelAB,gVSTel,
                       ncol=1,nrow=2,
                       heights=c(1.5,5),
                       labels=c('','F'),
                       align='v',
                       font.label = list(size=8),
                       hjust=0,vjust=1)

gAll <- ggarrange(gPanelABF,
                  gPanelSort,
                  ncol=2,nrow=1,
                  widths=c(2,1))

fScale <- 1
ggsave('Pratto_et_al_Figure7.png',plot = gAll, width=6.5*fScale,height=8*fScale,dpi = 300)
ggsave('Pratto_et_al_Figure7.pdf',plot = gAll, width=6.5*fScale,height=8*fScale)