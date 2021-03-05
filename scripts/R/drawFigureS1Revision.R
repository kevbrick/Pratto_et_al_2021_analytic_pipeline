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
library(dplyr)
library(magrittr)
library(ggupset)
library(cowplot)
library(patchwork)

theme7point()

## Function to plot mean signals

## HeatMaps
makeHMsForS1v2 <- function(set="WT45"){
  
  plotMNoverHM <- function(dfMn,dfHM,logScore=TRUE,facetSamples=TRUE,vHeight=c(2,3)){
    
    gMean <- ggplot(dfMn,
                    aes(x=pos+0.1,y=score,color=name)) +
      geom_vline(xintercept=0,lty='dashed',lwd=.1) +
      geom_smooth(lwd=.4,span=.2,se = FALSE) +
      scale_color_manual(values=c('dodgerblue2','firebrick','grey50')) +
      scale_x_continuous(breaks=c(-7,0,7)) +
      theme(legend.position=c(1,1),
            legend.justification=c(1,1),
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.key.size = unit(0.2,'cm'),
            axis.text.x=element_blank(),
            strip.background=element_blank(),
            strip.text = element_text(size=7),
            plot.margin = unit(c(0,0,0,0),'cm')) +
      ylab('Coverage') +
      xlab('') +
      coord_cartesian(expand=FALSE)
    
    if (logScore ){
      dfHM$dir <- 1 
      dfHM$dir[dfHM$score < 0] <- -1
      dfHM$score <- (dfHM$score ^ 2) * dfHM$dir
    }
    
    gHM <- ggplot(dfHM,
                  aes(x=pos,
                      y=hs,
                      fill=score)) + 
      geom_tile() + 
      theme(legend.position='none',
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.background=element_blank(),
            strip.text=element_blank(),
            panel.border=element_rect(size=.2,fill=NA),
            plot.margin = unit(c(0,0,0,0),'cm')) + 
      scale_x_continuous(breaks=c(-7,0,7)) +
      scale_fill_gradient2(low='red',mid='white',high='dodgerblue2',midpoint=0) + 
      ylab('Per locus coverage') + 
      xlab('Distance to center (Kb)') + 
      coord_cartesian(expand=FALSE)

    if (facetSamples){
      gMean <- gMean + facet_grid(~sample+dets+loci)
      gHM   <- gHM   + facet_wrap(label~sample+dets+loci,scales='free_y',nrow=2)
    }else{
      gHM   <- gHM   + facet_grid(label~.,scales='free_y') 
    }
    
    gX2 <- plot_grid(gMean,
                     gHM,
                     ncol=1,nrow=2,
                     align = 'v',
                     axis = 'lr')
    
    #gX2 <- ggarrange(apX2[[1]],apX2[[2]],ncol=1,heights=vHeight)
    
    return(list(gMean=gMean,
                gHM=gHM,
                gFig=gX2))
  }
  
  getHMplots <- function(matName,sample='s',dets='d',loci='name'){
    dfRet   <- plotFromDeepToolsMatrix(inMatrix=paste0(matName,'.fragsMatrix'),
                                       sortOrder = c('Watson','Crick','Both'),
                                       returnOrder = c('Crick','Watson','Both'),
                                       matrixType = 'referencepoint',
                                       samplesToExclude = 'Both',
                                       withHeatmap = TRUE,
                                       fSize = 10,winVE=100,winHZ=10)
    
    dfRet$figHM <- dfRet$figHM + geom_vline(xintercept=0,color='red',lwd=.1,lty='dashed') +
      scale_fill_gradient2(low='white',high='grey30') +
      scale_x_continuous(breaks=c(-5,0,5),labels=c(-5,0,5))
    
    dfRet$mData$sample <- sample
    dfRet$hmData$sample <- sample
    
    dfRet$mData$dets <- dets
    dfRet$hmData$dets <- dets
    
    dfRet$mData$loci <- loci
    dfRet$hmData$loci <- loci
    
    dfRet$hmData$score[dfRet$hmData$label == 'Watson'] <- dfRet$hmData$score[dfRet$hmData$label == 'Watson']*-1
    
    return(dfRet)
  }

  if (set == "WT12"){
    dfN1  <-  getHMplots('mm10_WT_Rep1.origins'        ,'WT1' , '-' , 'Ori')
    dfN2  <-  getHMplots('mm10_WT_Rep2.origins'        ,'WT2' , '-' , 'Ori')
    dfN3  <-  getHMplots('mm10_WT_Rep3.origins'        ,'WT3' , '-' , 'Ori')
    dfESC1 <-  getHMplots('mm10_ESC_Rep1.origins'      ,'ESC1', '-', 'Ori')
    dfESC2 <-  getHMplots('mm10_ESC_Rep2.origins'      ,'ESC2', '-', 'Ori')
    dfESR2 <-  getHMplots('mm10_RNAse_ESC_Rep2.origins','ESC2', '+RNAse', 'Ori')

    dfFacetMn <- rbind(dfN1$mData,
                       dfN2$mData,
                       dfN3$mData,
                       dfESC1$mData,
                       dfESC2$mData,
                       dfESR2$mData)

    dfFacetHM <- rbind(dfN1$hmData,
                       dfN2$hmData,
                       dfN3$hmData,
                       dfESC1$hmData,
                       dfESC2$hmData,
                       dfESR2$hmData)
  }
  
  if (set == 'WT45'){
    dfN4  <-  getHMplots('mm10_WT_Rep4.origins'       ,'WT4','-' ,'Ori')
    dfR4  <-  getHMplots('mm10_RNAh_WT_Rep4.origins'  ,'WT4','+RNAh' ,'Ori')
    dfN5  <-  getHMplots('mm10_WT_Rep5.origins'       ,'WT5','-' ,'Ori')
    dfR5  <-  getHMplots('mm10_RNAh_WT_Rep5.origins'  ,'WT5','+RNAh' ,'Ori')
    dfSp  <-  getHMplots('mm10_sperm_Rep1.origins'    ,'Sperm'  ,'-' ,'Ori')
  
    dfG4C4 <-  getHMplots('mm10_WT_Rep4.g4.Crick'     ,'WT4','-' ,'G4(C)')
    dfG4W4 <-  getHMplots('mm10_WT_Rep4.g4.Watson'    ,'WT4','-' ,'G4(W)')
    dfG4C5 <-  getHMplots('mm10_WT_Rep5.g4.Crick'     ,'WT5','-' ,'G4(C)')
    dfG4W5 <-  getHMplots('mm10_WT_Rep5.g4.Watson'    ,'WT5','-' ,'G4(W)')
    
    dfG4W4h <-  getHMplots('mm10_RNAh_WT_Rep4.g4.Crick'  ,'WT4','+RNAh' ,'G4(C)')
    dfG4C4h <-  getHMplots('mm10_RNAh_WT_Rep4.g4.Watson' ,'WT4','+RNAh' ,'G4(W)')
    dfG4WS  <-  getHMplots('mm10_sperm_Rep1.g4.Crick'  ,'Sperm','+RNAh' ,'G4(C)')
    dfG4CS  <-  getHMplots('mm10_sperm_Rep1.g4.Watson' ,'Sperm','+RNAh' ,'G4(W)')
    
    dfFacetMn <- rbind(dfN4$mData,
                       dfR4$mData,
                       dfN5$mData,
                       dfR5$mData,
                       dfSp$mData,
                       dfG4C4$mData,
                       dfG4W4$mData,
                       dfG4C4h$mData,
                       dfG4W4h$mData,                      
                       dfG4C5$mData,
                       dfG4W5$mData,
                       dfG4CS$mData,
                       dfG4WS$mData)
    
    dfFacetHM <- rbind(dfN4$hmData,
                       dfR4$hmData,
                       dfN5$hmData,
                       dfR5$hmData,
                       dfSp$hmData,
                       dfG4C4$hmData,
                       dfG4W4$hmData,
                       dfG4C4h$hmData,
                       dfG4W4h$hmData,                      
                       dfG4C5$hmData,
                       dfG4W5$hmData,
                       dfG4CS$hmData,
                       dfG4WS$hmData)
  }
  
  #sOrder <- c('WT Rep4','WT Rep5','WT Rep1','WT Rep2','WT Rep3','Sperm','ESC Rep1','ESC Rep2')
  sOrder <- c('WT4','WT5','WT1','WT2','WT3','Sperm','ESC1','ESC2')
  dfFacetMn$sample <- factor(dfFacetMn$sample,levels=sOrder)
  dfFacetHM$sample <- factor(dfFacetHM$sample,levels=sOrder)
  
  detsOrder <- c('-','-RNAh','+RNAh','-RNAse','+RNAse')
  dfFacetMn$dets <- factor(dfFacetMn$dets,levels=detsOrder)
  dfFacetHM$dets <- factor(dfFacetHM$dets,levels=detsOrder)
  
  lociOrder <- c('Ori','G4(C)','G4(W)')
  dfFacetMn$loci <- factor(dfFacetMn$loci,levels=lociOrder)
  dfFacetHM$loci <- factor(dfFacetHM$loci,levels=lociOrder)
  
  dfPlots <-  plotMNoverHM(dfFacetMn,dfFacetHM,logScore = FALSE)
  
  return(dfPlots)
}

getOriginOverlaps <- function(dfOL){
  ori_samples <- names(dfOL)[grepl('^ori_',names(dfOL)) & !grepl('(RNAse|sperm|RNAh|Rep4|Rep5)',names(dfOL))]
  ori_names   <- gsub('ori_','',names(dfOL)[grepl('^ori_',names(dfOL)) & !grepl('(RNAse|sperm|RNAh|Rep4|Rep5)',names(dfOL))])

  dfAnother <- data.frame(name=ori_names,
                          count=0,
                          shared=0,
                          unique=0,
                          pc=0,
                          label="")
  
  for (s1 in ori_samples){
    name1 <- gsub('ori_','',s1)
    dfAnother$shared[dfAnother$name == name1] <- sum(dfOL[[s1]]>0 & dfOL$oriSum>1)
    dfAnother$unique[dfAnother$name == name1] <- sum(dfOL[[s1]]>0 & dfOL$oriSum==1)
    dfAnother$count[dfAnother$name == name1]  <- sum(dfOL[[s1]]>0)
    dfAnother$pc[dfAnother$name == name1]     <- dfAnother$shared[dfAnother$name == name1] / dfAnother$count[dfAnother$name == name1] * 100 
    dfAnother$label[dfAnother$name == name1]  <- round(dfAnother$pc[dfAnother$name == name1],0)
  }
  
  dfAnother <- dfAnother %>% mutate(namelbl=paste0(gsub("_"," ",name),
                                                   "\n(N =",format(count,big.mark = ","),
                                                   ")"))
  
  mAnother <- reshape2:::melt.data.frame(dfAnother,
                                   id.vars = 'namelbl', 
                                   measure.vars = c('unique','shared'),
                                   variable.name = 'type',
                                   value.name='count')

  gOL <- ggplot(mAnother) + 
    geom_bar(aes(x=count/1000,
                 y=namelbl,
                 fill=type),
             position='stack',stat='identity',
             lwd=.3,color='black') + 
    geom_text(data=dfAnother,
              aes(y=namelbl,
                  label=paste0(round(pc,0),'%')),
              x=1,
              size=7*5/14,
              hjust=0) + ylab('') + 
    xlab(bquote('Ori-SSDS peaks (x'*10^3*')'))+ 
    theme(legend.position=c(1,.4),
          legend.key.size=unit(0.2,'cm'),
          legend.justification=c(1,.5),
          legend.direction='vertical',
          legend.title=element_blank(),
          legend.text=element_text(size=7),
          legend.background=element_blank(),
          axis.text.y=element_text(hjust=0.5)) + 
    scale_fill_manual(values=c('pink','grey80')) + 
    coord_cartesian(xlim=c(0,19),ylim=c(.4,5.6),expand=FALSE)
  
  return(gOL)
}

### MAIN FUNCTION
lOri <- getOriginsDFMM10('mm10_OriSSDS.Origins.tab')

allOriDF <- lOri$ori

########## SET 1: ALL ORIS in ANY WT SAMPLE : REQUIRE W/C ASYMMETRY OR BE >7Kb
oriUnionDF <- parseOriginsMM10(dfOri = allOriDF,oType='union')

########## SET 2: FINAL "HIGH CONFIDENCE" ORIS; IN >1 SAMPLE : REQUIRE W/C ASYMMETRY OR BE >7Kb
oriDF      <- parseOriginsMM10(allOriDF,oType='hiconf')

########## SCATTERPLOTS
gSc        <- plotScattersFigMM10(oriUnionDF$ori,oriDF$ori)

########## ESCs
dfESC <- parseOriginsTestis_V_ESCs(allOriDF)

################## UpSet Plot
dfUpset <- dfESC %>% 
  mutate("WT Rep 1"=ori_WT_Rep1>0,
         "WT Rep 2"=ori_WT_Rep2>0,
         "WT Rep 3"=ori_WT_Rep3>0,
         "ESC Rep 1"=ori_ESC_Rep1>0,
         "ESC Rep 2"=ori_ESC_Rep2>0,
         name=paste(cs,from,to,sep=":")) %>% 
  select(name,
         ESC1="ESC Rep 1",
         ESC2="ESC Rep 2",
         WT1="WT Rep 1",
         WT2="WT Rep 2",
         WT3="WT Rep 3") %>% 
  reshape2:::melt(id='name',value.name='ok', variable.name = 'sample') %>% 
  dplyr:::filter(ok) %>% 
  select(name,sample) 
  
dfUpsetPlot <- dfUpset %>% 
  group_by(name) %>% 
  summarize(samples = list(sample),slst=paste(list(sample))) %>% 
  group_by(slst) %>% add_count(name="grpCount") %>%
  dplyr:::filter(grpCount > 50) %>%
  select(name,samples)

gUpsetOnly  <- ggplot(dfUpsetPlot,
                  aes(x = samples)) +
  geom_bar() +
  scale_x_upset() + xlab('') + ylab('Origins (thousands)') + 
  geom_text(stat='count', 
            aes(label=format(after_stat(count),big.mark = ",")), 
            hjust=-0.1,vjust=-1,size=7*5/14,angle=45) + 
  geom_blank(aes(y=10000)) + 
  scale_y_continuous(breaks=c(0,5000,10000), labels=c(0,5,10)) + 
  theme_combmatrix(combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.line.size = .5)

gOverlaps <- getOriginOverlaps(dfESC)
gUpset    <- gUpsetOnly + inset_element(ggarrange(gOverlaps,labels='F',font.label = list(size=8,fontface='bold'),
                                                  hjust=0,vjust=1),0.2,0.3,1,1)

## Scatterplots V ESC Origins
gESC  <- plotScatters_Testis_V_ESCs(dfESC)

gJKL <- ggarrange(gESC$gESC2WT,
                  gESC$gEEu,
                  gESC$gEEi,
                  ncol=1,nrow=3,
                  labels=c('J','K','L'),
                  font.label = list(size=8,fontface='bold'),
                  hjust=0,vjust=1)

gGHIJKL <- ggarrange(gSc$gV1,
                     gSc$gV2,
                     ggplot()+theme_void(), 
                     gESC$gVWTc,
                     gJKL,
                     ncol=5,
                     nrow=1,
                     align='h',
                    widths=c(3.2,3,.05,3,3),
                     labels=c('G','H','','I',''),
                     font.label = list(size=8,face='bold',hjust=0,vjust=1))

gEtoL <- ggarrange(gUpset,
                   gGHIJKL,
                   nrow=1,
                   ncol=2,
                   widths=c(9,12),
                   labels=c('E',''),
                   font.label = list(size=8,face='bold',hjust=0,vjust=1))

## Get Slice
mySlice <- 'chr2_161Mb_figureS1'
sliceWT <- drawOriSlice(slice = mySlice,
                        slicefolder=slicedir,
                        pattern2omit = 'ep[123]',
                        winType='w2000s100',
                        noATAC = TRUE,
                        noCGI=TRUE,
                        noG4=TRUE,
                        filled=FALSE,
                        annotScale=0.35)

sliceES2 <- drawOriSlice(slice = mySlice,
                        samples2use = c('WT_Rep1','WT_Rep2','ESC_Rep1','ESC_Rep2','RNAse_ESC_Rep2'),
                        winType='w2000s100',
                        noATAC = TRUE,
                        noCGI=TRUE,
                        noG4=TRUE,
                        filled=FALSE,
                        annotScale=0.35)

gSlicesAligned <- align_plots(sliceES2$gCover+ theme(plot.margin=unit(c(0,0,0,0),'cm')),
                              sliceWT$gCover+ theme(plot.margin=unit(c(0,0,0,0),'cm')),
                              sliceWT$gAnnot+ theme(plot.margin=unit(c(0,0,0,0),'cm')),
                              align='v',axis = 'lr')

gSlicesAll <- ggarrange(gSlicesAligned [[1]],
                     gSlicesAligned [[2]],
                     gSlicesAligned [[3]],
                     nrow=3,heights=c(4,4,2),
                     align='v')

gOK <-  rasterGrob(readPNG('OKSeq.HM.png'), interpolate=TRUE)
                   
gSlices <- ggarrange(gSlicesAll,
                     gOK,
                     ncol=2,widths=c(2,1),
                     labels=c('','D'),
                     font.label = list(size=8,face='bold'),
                     hjust=0,vjust=1)

# ## Heatmaps and mean cover
gAveragePlots12 <- makeHMsForS1v2("WT12")
gAveragePlots45 <- makeHMsForS1v2("WT45")

gAveragePlots <- ggarrange(gAveragePlots12$gFig,
                           gAveragePlots45$gFig,
                           ncol=2,nrow=1,
                           widths=c(2,3),
                           labels=c('B','C'),
                           font.label = list(size=8,face='bold'),
                           hjust=0,vjust=1)

gAll <- ggarrange(gSlices,
                  gAveragePlots,
                  gEtoL,
                  nrow=3,
                  heights=c(5,4,5),
                  labels=c('A','',''),
                  font.label = list(size=8,face='bold'),
                  hjust=0,vjust=1)

ggsave(paste0('Pratto_Et_Al_FigS1_',mySlice,'.png'),gAll,width=7.5,height=10)
ggsave(paste0('Pratto_Et_Al_FigS1_',mySlice,'.pdf'),gAll,width=7.5,height=10)

