imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(data.table)
library(ggrepel)

theme7point()

drawSNScomparisonFigs <- function(){
  
  dfSNS <- read.table('snsClusters.tab')
  names(dfSNS) <- c('src','sz')
  
  ## Get counts by cluster size
  dfSNSPlot <- dfSNS     %>% group_by(src,sz)   %>% 
    count(name='count') %>% mutate(nOri=sz*count, 
                                   ref=gsub("_"," ",src),
                                   src=ifelse(grepl('Cayrou',src),'Mouse','Human'))
  
  ## Get origins vs clusters counts for barplot
  dfSNSBars <- dfSNSPlot %>% group_by(src,sz>0) %>% 
    add_tally(nOri,name='SNS-Seq peaks') %>%
    group_by(src) %>% 
    add_tally(count, name='Peak clusters') %>% 
    dplyr:::filter(sz==1) %>%
    select(src,ref,'SNS-Seq peaks','Peak clusters') %>% 
    reshape2:::melt(value.name = "N", variable.name = "type")
  
  ## Fix order
  dfSNSBars$ref  <- factor(dfSNSBars$ref,levels=c('Long 2020','Cayrou 2015'))
  dfSNSBars$src  <- factor(dfSNSBars$src,levels=c('Human','Mouse'))
  dfSNSBars$type <- factor(dfSNSBars$type,levels=c('Peak clusters','SNS-Seq peaks'))
  
  gSNSBarsh<- ggplot(dfSNSBars,aes(x=src,y=N/1000,fill=type)) +
    geom_bar(color='black',lwd=.3,stat='identity',position=position_dodge()) +
    xlab('') + ylab(bquote('Origin count (x'*10^3*')')) +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background=element_blank(),
          legend.key.size=unit(0.3,'cm'),
          legend.text=element_text(size=7)) +
    geom_text(aes(label=paste0("N=",round(N/1000)," "),y=N/1000), 
              position=position_dodge(width=1),
              hjust=1,
              size=7*5/14) +
    coord_flip() +
    scale_fill_manual("",values = c('firebrick','pink'),guide = guide_legend(reverse=TRUE))

  gSNSBars<- ggplot(dfSNSBars,aes(x=src,y=N/1000,fill=type)) +
    geom_bar(color='black',lwd=.3,width=.9, stat='identity',position=position_dodge()) +
    xlab('') + ylab(bquote('Origin count (x'*10^3*')')) +
    theme(legend.position=c(0,1),
          legend.justification=c(0,1),
          legend.background=element_blank(),
          legend.key.size=unit(0.2,'cm'),
          legend.text=element_text(size=7),
          axis.text.y=element_text(angle=90,hjust=.5),
          axis.ticks.y=element_blank()) +
    geom_text(aes(label=paste0("N=",round(N/1000)," "),y=(N/1000)-1), 
              position=position_dodge(width=.9),
              vjust=1,
              hjust=0.5,
              size=7*5/14) +
    coord_cartesian(xlim=c(.54,2.46),expand=FALSE) +
    scale_fill_manual("",values = c('firebrick','pink'),guide = guide_legend(reverse=TRUE))
  
  ## Pies
  dfSNSPie <- dfSNSPlot %>% group_by(src) %>% 
    mutate(type=ifelse(sz>1,'In cluster','Isolated')) %>%
    add_tally(nOri, name='total') %>% 
    group_by(src,total,type) %>% 
    tally(nOri,name='N') %>% 
    mutate(pc=N/total*100,label=paste0(round(pc,0),'%'),
           labely=ifelse(pc > 50,95,5),
           yjust=ifelse(pc > 50,1,0))
  
  dfSNSPie$type <- factor(dfSNSPie$type,levels=c('In cluster','Isolated'))
  
  gSNSPies <- ggplot(dfSNSPie,aes(x=1,y=pc,fill=type)) +
    geom_bar(stat='identity',color='black',lwd=.3) +
    coord_polar("y",start=0) +
    geom_label(aes(y=pc-10,label = label), 
               size=7*5/14,color='black',
               fill='white') +
    theme(axis.text=element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position='bottom',
          legend.title=element_text(size=10),
          legend.spacing.x=unit(0.25,'cm'),
          legend.key.size = unit(0.25,'cm'),
          strip.background=element_blank(),
          strip.text = element_text(size=7),
          plot.margin = margin(0,0,0,0,'cm')) +
    facet_grid(src~.,switch='y') +
    xlab('') + ylab('') +
    scale_fill_manual("",values = c('grey20','grey50'))

  gSNSStacks <- ggplot(dfSNSPie,aes(x=src,y=pc,fill=type)) +
    geom_bar(stat='identity',width=.8,color='black',lwd=.3) +
    coord_flip(expand=FALSE,xlim=c(.5,2.5)) + 
    geom_label(aes(y=labely,hjust=yjust,label = label), 
               size=7*5/14,color='black',
               fill='white') +
    theme(axis.ticks = element_blank(),
          legend.position='top',
          legend.title=element_text(size=10),
          legend.spacing.x=unit(0.25,'cm'),
          legend.key.size = unit(0.25,'cm'),
          legend.text = element_text(size=7),
          strip.background=element_blank(),
          strip.text = element_text(size=7),
          plot.margin = margin(0,0,0,0,'cm')) +
    xlab('') + ylab('Origins (%)') + 
    scale_fill_manual("",values = c('grey20','grey50'))
  
  gSNSStacksHZ <- ggplot(dfSNSPie,aes(x=src,y=pc,fill=type)) +
    geom_bar(stat='identity',width=0.9,color='black',lwd=.3) +
    geom_text(aes(y=labely,vjust=yjust,label = label), 
               size=7*5/14,color='white') +
    theme(axis.ticks = element_blank(),
          axis.text.x=element_text(angle=25,hjust=1),
          legend.position='top',
          legend.direction='vertical',
          legend.title=element_text(size=10),
          legend.spacing.x=unit(0.1,'cm'),
          legend.spacing.y=unit(0.2,'cm'),
          legend.key.size = unit(0.2,'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          legend.text = element_text(size=7),
          strip.background=element_blank(),
          strip.text = element_text(size=7),
          plot.margin = margin(0,0,0,0,'cm')) +
    xlab('') + ylab('Origins (%)') + 
    coord_cartesian(xlim=c(0.5,2.5),expand=FALSE) +
    scale_y_continuous(breaks=c(0,100)) + 
    scale_fill_manual("",values = c('grey20','grey50'))
  
  ## Import example image
  mypng <- readPNG(source = 'Pratto_et_al_Figure2A_SNS_v_OriSSDS_one_origin.png')
  gPic <- rasterGrob(mypng, interpolate=TRUE)
  
  return(list(gExample=gPic,
              gPies=gSNSPies,
              gPieStack=gSNSStacks,
              gPieStackHZ=gSNSStacksHZ,
              gBars=gSNSBars))
}

makeHMs <- function(){

  plotMNoverHM <- function(dfMn,dfHM,logScore=TRUE,facetSamples=TRUE,vHeight=c(1,3),noLeg=FALSE){
    
    gMean <- ggplot(dfMn,
                    aes(x=pos+0.2,y=score,color=name)) +
      geom_vline(xintercept=0,lty='dashed',lwd=.1) +
      geom_smooth(lwd=.4,span=.2,se = FALSE) +
      scale_color_manual(values=c('firebrick','dodgerblue2','grey50')) +
      scale_x_continuous(breaks=c(-3,0,3)) +
      theme(legend.position=c(1,1),
            legend.justification=c(1,1),
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.key.size = unit(0.2,'cm'),
            axis.text.x=element_blank(),
            strip.background=element_blank(),
            strip.text = element_text(size=7),
            plot.margin = margin(0,0,0,0,'cm')) +
      ylab('Coverage') +
      xlab('') +
      coord_cartesian(expand=FALSE)
    
    if (logScore ){
      dfHM$dir <- 1 
      dfHM$dir[dfHM$score < 0] <- -1
      dfHM$score <- (dfHM$score ^ 2) * dfHM$dir
    }
    
    gHM <- ggplot(dfHM,
                  aes(x=pos+.2,
                      y=hs,
                      fill=score)) + 
      geom_tile() + 
      theme(legend.position='none',
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.background=element_blank(),
            strip.text=element_blank(),
            panel.border=element_rect(size=.2,fill=NA),
            plot.margin = margin(0,0,0,0,'cm')) + 
      scale_x_continuous(breaks=c(-3,0,3)) +
      scale_fill_gradient2(low='red',mid='white',high='dodgerblue2',midpoint=0) + 
      ylab('Per locus coverage') + 
      xlab('Distance to\nCpG peak (Kb)') + 
      coord_cartesian(expand=FALSE)
    
    if (facetSamples){
      gMean <- gMean + facet_wrap(~sample,nrow=1) 
      gHM   <- gHM   + facet_wrap(label~sample,scales='free_y',nrow=2)
    }else{
      gHM   <- gHM   + facet_grid(label~.,scales='free_y') 
    }
    
    if (noLeg){
      gMean <- gMean + theme(legend.position='none')
      gHM   <- gHM   + theme(legend.position='none')
    }
    
    apX2 <- align_plots(gMean,
                        gHM,
                        align = 'v',axis = 'lr')
    
    gX2 <- ggarrange(apX2[[1]],apX2[[2]],ncol=1,heights=vHeight)
    
    return(list(gMean=gMean,
                gHM=gHM,
                gFig=gX2))
  }
  
  getHMplots <- function(sample='Sperm',dets='NA',loci='Origins'){
    
    sampleExclude <- ifelse(sample=='Sperm','Testis','Sperm')
    
    dfRet <- plotFromDeepToolsMatrix(inMatrix=paste0('CpG_atOriBroad.peaks.Matrix'),
                                 sortOrder = c('Crick (Testis)','Watson (Testis)',
                                               'Crick (Sperm)','Watson (Sperm)'),
                                 matrixType = 'referencepoint',
                                 samplesToExclude = paste0(c('Crick','Watson'),paste0(' (',sampleExclude,')')),
                                 normByRow = TRUE,
                                 winVE = 8,
                                 withHeatmap = TRUE)
    
    dfRet$figHM <- dfRet$figHM + geom_vline(xintercept=0,color='red',lwd=.1,lty='dashed') +
      scale_fill_gradient2(low='white',high='grey30') +
      scale_x_continuous(breaks=c(-5,0,5),labels=c(-5,0,5))
    
    dfRet$mData$sample <- sample
    dfRet$hmData$sample <- sample
    
    dfRet$mData$dets <- dets
    dfRet$hmData$dets <- dets
    
    dfRet$mData$loci <- loci
    dfRet$hmData$loci <- loci

    dfRet$mData$name <- factor(dfRet$mData$name, levels=c('Watson','Crick', paste0(c('Crick','Watson'),paste0(' (',sample,')'))))
    dfRet$mData$name[grepl('Watson',dfRet$mData$name)] <- 'Watson'
    dfRet$mData$name[grepl('Crick',dfRet$mData$name)] <- 'Crick'
    
    dfRet$hmData$label <- factor(dfRet$hmData$label, levels=c('Watson','Crick', paste0(c('Crick','Watson'),paste0(' (',sample,')'))))
    dfRet$hmData$label[grepl('Watson',dfRet$hmData$label)] <- 'Watson'
    dfRet$hmData$label[grepl('Crick',dfRet$hmData$label)] <- 'Crick'
      
    dfRet$hmData$score[dfRet$hmData$label == 'Watson'] <- dfRet$hmData$score[dfRet$hmData$label == 'Watson']*-1
    
    return(dfRet)
  }
  
  dfSperm  <-  getHMplots('Sperm',  'NA', 'Origins')
  dfWT     <-  getHMplots('Testis', 'NA', 'Origins')

  dfFacetMn <- rbind(dfSperm$mData,
                     dfWT$mData)
  
  dfFacetHM <- rbind(dfSperm$hmData,
                     dfWT$hmData)
  
  sOrder <- c('Testis','Sperm')
  dfFacetMn$sample <- factor(dfFacetMn$sample,levels=sOrder)
  dfFacetHM$sample <- factor(dfFacetHM$sample,levels=sOrder)
  
  dfPlots <-  plotMNoverHM(dfFacetMn,dfFacetHM,logScore = FALSE)
  
  return(list(gSperm  = plotMNoverHM(dfSperm$mData,dfSperm$hmData,logScore = FALSE, facetSamples = FALSE, noLeg=TRUE),
              gTestis = plotMNoverHM(dfWT$mData,dfWT$hmData,logScore = FALSE, facetSamples = FALSE, noLeg=TRUE),
              gX2     = plotMNoverHM(dfFacetMn,dfFacetHM,logScore = FALSE)))
}

drawOriSSDSvCpGvClusters <- function(){
  
  ### Plot CpG peaks v origin center @ narrow origins
  dfAllOriBed  <- read.table('ori7k.bed')
  dfCpGOriBed  <- read.table('oneCpg.bed')
  
  nTot         <- dim(dfAllOriBed)[1]
  nCpG         <- dim(dfCpGOriBed)[1]
  
  dfCpGDist    <- read.table('oriVCpG.tab') %>% 
    mutate(distance=round(V1/20)*20/1000) %>% 
    group_by(distance) %>% count(name='N') %>%
    mutate(freq=N/dim(dfCpGOriBed)[1] * 100)

  gCpG_v_Ori <- ggplot(dfCpGDist,aes(x=distance,y=freq)) +
    geom_vline(lwd=.1,lty='dashed',xintercept=c(-3,-2,-1,1,2,3)) +
    geom_vline(lwd=.1,xintercept=0) +
    geom_point(size=0.4) +
    geom_smooth(span=.2,lwd=.2,fill=NA) +
    xlab('Distance to origin center (Kb)') +
    ylab('CpG peaks (%)') +
    theme(plot.title = element_text(hjust=0.5)) + 
    coord_cartesian(expand=FALSE,ylim=c(-0.2,3))
  ##ggtitle(annotText) +  
  
  ### Plot CpGs Vs W/C Heatmaps
  gHeatMaps <- makeHMs()
  
  ## Origin width plot
  dfOriWidth <- read.table('hiconf_origins.mm10.bedgraph') %>% 
    mutate(cs=V1,from=V2,to=V3,nOri=V4, width=ceiling((to-from)/1000)) %>%
    group_by(width) %>% count(name='N') %>%
    ungroup() %>%
    mutate(pc=cumsum(N)/max(cumsum(N))*100)
  
  gCumSum <- ggplot(dfOriWidth,aes(x=width,y=pc)) +
    geom_point(size=.4) +
    geom_line(lwd=.2) +
    scale_x_log10() +
    xlab('Origin width (Kb)') +
    ylab('Cumulative\norigins (%)') +
    annotation_logticks(sides='b',
                        short=unit('0.1','cm'),
                        mid=unit('0.1','cm'),
                        long=unit('0.2','cm'),
                        lwd=.2,size=.2) + 
    coord_cartesian(xlim=c(1,110),expand=FALSE)
  
  ## Origin width plot
  dfOriWidth <- read.table('CpGpeaksAtOrigins.txt') %>% 
    mutate(cs=V1, from=V2, to=V3, nCpG=V4, width=to-from,
           type=ifelse(width < 5000,'0 - 5',
                       ifelse(width < 10000,'5 - 10',
                              ifelse(width < 15000,'10 - 15', 
                                     ifelse(width < 20000,'15 - 20', '> 20')))),
           CpG_count=ifelse(nCpG>3,'4+',nCpG)) %>%
    group_by(type) %>% add_count(name='total') %>%
    group_by(type,CpG_count,total) %>% count(name='N') %>%
    mutate(pc=N/total*100)
  
  dfOriWidth$type <- factor(dfOriWidth$type, levels=c('0 - 5','5 - 10','10 - 15','15 - 20', '> 20'))
  
  barCols   <- c('0'  = "#E69F00", 
                 '1'  = "#56B4E9", 
                 '2'  = "#D55E00", 
                 '3'  = "#444444",
                 '4+' = "#cccccc")
  
  gCpgAtOriBar <- ggplot(dfOriWidth,aes(x=type,y=pc,fill=factor(CpG_count))) +
    geom_bar(stat='identity',color='black',lwd=.3) +
    coord_flip() + 
    xlab('Origin width (Kb)') + 
    ylab('Origins with N CpG peaks (%)') +
    theme(legend.key.size = unit(0.2,'cm'),
          legend.position = 'top') + 
    scale_fill_manual('# CpG peaks',values=barCols)
  
  #### Regions of interest
  gSlice_HoxA    <- drawOriClusterSlice(slice='chr6_52Mb_figure2', samples2use = c('sperm_Rep1','WT_Rep1','WT_Rep2','WT_Rep3'),useOrder=TRUE)
  gSlice_Lots    <- drawOriClusterSlice(slice='chr2_74Mb_figure2', samples2use = c('sperm_Rep1','WT_Rep1','WT_Rep2','WT_Rep3'),useOrder=TRUE)
  #gSlice_1_0     <- drawOriClusterSlice(slice='chr15_84Mb_1_0_twoWide',samples2use = c('sperm_Rep1','WT_Rep1','WT_Rep2','WT_Rep3'),useOrder=TRUE)
  #gSlice_Complex <- drawOriClusterSlice(slice='chr15_74Mb_complex',    samples2use = c('sperm_Rep1','WT_Rep1','WT_Rep2','WT_Rep3'),useOrder=TRUE)
  
  return(list(gSlice1=gSlice_HoxA,
              gSlice2=gSlice_Lots,
              gOriWidth=gCumSum,
              gCpGDist=gCpG_v_Ori,
              gCpGAtOri=gCpgAtOriBar,
              gCpGHeatmap=gHeatMaps$gTestis))
}

lstSNS <- drawSNScomparisonFigs()
lstCpG <- drawOriSSDSvCpGvClusters()

### Draw figure
apX2 <- align_plots(lstCpG$gCpGHeatmap$gMean + scale_y_continuous(breaks=c(5,10,15),labels=c(5,'',15)),
                    lstCpG$gCpGHeatmap$gHM,
                    align = 'v',axis = 'lr')

gHM <- ggarrange(apX2[[1]],apX2[[2]],ncol=1,heights=c(1,2.5))

gBC <- ggarrange(lstSNS$gBars,
                 lstCpG$gOriWidth,
                 ggplot()+theme_void(),
                 gHM,
                 ncol=4,nrow=1,
                 widths=c(1,1,0.1,.7),
                 labels = c('C','D','F',''),
                 font.label = list(size=8,fontface='bold'),
                 hjust=0,vjust=1)

gDF <- ggarrange(lstCpG$gCpGDist,
                 lstCpG$gCpGAtOri,
                 ncol=2,nrow=1,
                 widths=c(1,1),
                 labels = c('E','G'),
                 font.label = list(size=8,fontface='bold'),
                 hjust=0,vjust=1)

gBCDEF <- ggarrange(gBC,
                   gDF,
                   ncol=1,nrow=2,
                   heights=c(1,1))


gTop <- ggarrange(lstSNS$gExample,gBCDEF,
                  ncol=2,nrow=1,
                  widths=c(5,12))

gBott <- ggarrange(lstCpG$gSlice1,
                   lstCpG$gSlice2,
                   ncol=2,nrow=1,
                   widths=c(1,1),
                   labels = c('H','I'),
                   font.label = list(size=8,fontface='bold'),
                   hjust=0,vjust=1)

gAll <- ggarrange(gTop,gBott,
                  ncol=1,nrow=2,
                  heights=c(3,3))

nScale <- 1
ggsave('Pratto_et_al_Figure2.png',plot = gAll, height=6*nScale,width=6.5, dpi = 400)
ggsave('Pratto_et_al_Figure2.pdf',plot = gAll, height=6*nScale,width=6.5)

nScale <- 1
ggsave('Pratto_et_al_Figure2_top.pdf'   ,plot = gTop, height=3*nScale,width=6.5)
ggsave('Pratto_et_al_Figure2_bottom.pdf',plot = gBott, height=3*nScale,width=6.5)

