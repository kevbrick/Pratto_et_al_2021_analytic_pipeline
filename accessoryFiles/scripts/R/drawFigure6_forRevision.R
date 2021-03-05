imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')

library(data.table)

theme7point()

library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)
library(ggpmisc)
library(ggpubr)

## FUNCTIONS
plot_replication_speed <- function(dfX,pc,showLegend=FALSE,mainFig=FALSE){
  
  dfX <- dfX[dfX$pc == pc,]
  ylbl <- paste0("Time to replicate ",pc,"% (RT-Sim cycles)")
  xLoVal <- 0
  
  if (mainFig){
    dfX$csType <- 'others'
    dfX$csType[dfX$cs == 11] <- 'chr11'
    dfX$csType[dfX$cs == 13] <- 'chr13'
    if (pc == 25){
      xLoVal <- 0
      xHiVal <- 250
    }
    
    if (pc == 50){
      xLoVal <- 100
      xHiVal <- 350
    }
    
    if (pc == 75){
      xLoVal <- 200
      xHiVal <- 450
    }
    
    legendX <- 0
    legendY <- 1
  }else{
    dfX$csType <- 'others'
    #dfX$csType[dfX$cs == 11] <- 'chr11: 122 Mb: 1,033 origins'
    dfX$csType[dfX$cs == 11] <- 'chr11: 8.5 ori/Mb'
    #dfX$csType[dfX$cs == 13] <- 'chr13: 120 Mb:   327 origins'
    dfX$csType[dfX$cs == 13] <- 'chr13: 2.7 ori/Mb'
    
    xLoVal <- 0
    xHiVal <- 450
    
    legendX <- 1
    legendY <- 1
  }
  
  g <- ggplot(dfX,aes(x=cs,y=time,group=cs,fill=csType)) + 
    geom_boxplot(outlier.alpha=0,lwd=.1) + 
    coord_flip(ylim=c(xLoVal,xHiVal)) +
    scale_fill_manual("",values=c('orange','red','grey65')) + 
    xlab('Chromosome') + 
    ylab(ylbl) + 
    scale_x_reverse(breaks=c(1,10,19)) + 
    scale_y_continuous(breaks=seq(0,450,by=100))
  
  if (showLegend){
    g <- g + theme(legend.position=c(legendX,legendY),
                   legend.justification = c(legendX,legendY),
                   legend.background=element_blank(),
                   legend.key.size = unit(0.3,'cm'),
                   legend.text = element_text(size = 7))
  }else{
    g <- g + theme(legend.position='none')
  }
  
  return(g)
}

getReplicationSpeed <- function(modelFile = 'MeiS_ALL_S3_2to4C_STRA8_DMRT1_mm10.Rdata', altDF=NULL){
  
  if (is.null(altDF)){
    #### Replication speed per chrom: get from RT-Sim model
    load(modelFile)
    
    ## This is a df that details the time at which each segment of each chromosome replicates
    ## n x m : n = genome pos; m = cell; value = replicaiton cycle
    repPerCell <- modelData$mod$repPerCell
    repCS      <- modelData$mod$model$whatCS
    
    ## 0 == not replicated ... set to 501 (max)
    repPerCell[repPerCell==0] <- 501
  }else{
    repPerCell <- altDF$repPerCell
    repCS      <- altDF$repCS
  }
  
  chroms <- c(1:19)
  
  ## Just a flag for the first pass
  pos = 1
  
  ## Loop through chroms:
  for (cs in chroms){
    print(paste0("chr",cs))
    csPos  <- which(repCS == paste0('chr',cs))
    csSz   <- length(csPos)
    thisCS <- repPerCell[csPos,]
    
    ## Loop through each cell (n = 100)
    for (i in 1:100){
      
      ## loop through time & note percent replicated @ each time point
      prevRepPC <- 0
      for (t in 1:501){
        
        ## How many positions were replicated at or before time t (as a % of chromosome)
        repPC <- round(sum(thisCS[,i] <= t)/csSz*100)
        
        ## This assures that the only the earliest time point is used for each %
        if (repPC > prevRepPC){
          
          dfThis <- data.frame(cs=cs,cell=i,pc=repPC,time=t)
          if (pos == 1){
            dfCS <- dfThis
          }else{
            dfCS <- rbind(dfCS,dfThis)
          }
          #dfCS$replicated[pos] <- t
          pos = pos+1
        }
        
        prevRepPC <- repPC
      }
    }
  }
  
  ## No need for every single time-point
  dfReduced <- dfCS[dfCS$pc %in% c(10,25,50,75,90),]
  
  theme7point()
  
  ## Set names for graph
  dfReduced$csType <- 'others'
  #  dfReduced$csType[dfReduced$cs == 11] <- 'chr11: 122 Mb: 1,033 origins'
  #  dfReduced$csType[dfReduced$cs == 13] <- 'chr13: 120 Mb:   327 origins'
  dfReduced$csType[dfReduced$cs == 11] <- 'chr11: 8.5 ori/Mb'
  dfReduced$csType[dfReduced$cs == 13] <- 'chr13: 2.7 ori/Mb'
  
  dfReduced$lbl <- paste0('Cycles to replicate ',dfReduced$pc,'% of chromosome')
  
  save(dfReduced,file='replicationSpeedByCS.Rdata')
  #load(file = 'replicationSpeedByCS.Rdata')
  
  gMain <- plot_replication_speed(dfReduced,75,showLegend = TRUE, mainFig = TRUE)
  
  g25     <- plot_replication_speed(dfReduced,25,showLegend = TRUE)
  g50     <- plot_replication_speed(dfReduced,50)
  g75     <- plot_replication_speed(dfReduced,75)
  
  gX3 <- ggarrange(g25,g50,g75,ncol=1,nrow=3)
  
  return(list(gMain = gMain,
              gX3   = gX3,
              g25   = g25,
              g50   = g50,
              g75   = g75,
              df    = dfReduced))
}

plotObsVsPred <- function(dfPlot, res="COperMb", preds=preds, name='xx'){
  
  f <- as.formula( paste(res, paste(preds, collapse = " + "), sep = " ~ "))
  
  fit <- eval(bquote(lm(.(f), data = dfPlot )) )
  summ <- summary(fit)
  
  dfPlot$residual <- resid(fit)
  dfPlot$prediction <- predict(fit, newdata = dfPlot)
  
  dfPlot$selectCS   <- '' 
  dfPlot$selectCS[dfPlot$chr == 'chr11']   <- 'chr11'
  dfPlot$selectCS[dfPlot$chr == 'chr13']   <- 'chr13'
  
  dfPlot$selectCSv   <- 0 
  dfPlot$selectCSv[dfPlot$chr == 'chr11']   <- -1.5
  dfPlot$selectCSv[dfPlot$chr == 'chr13']   <- 2
  
  #r2Label <- bquote("   "*R^2*" = "*.(round(summ$r.squared,2)))
  
  g <- ggplot(dfPlot, aes(x = prediction, y = COperMb)) +
    geom_smooth(method = "lm", se=T,alpha=0.25, lwd=.4, colour='grey') +
    geom_point(aes(fill = selectCS, color = selectCS), size=0.8) +
    annotate(geom='text',label=paste0(" ",name),x=Inf,y=-Inf,size=7*5/14, hjust=1, vjust=-1) + 
    annotate(geom='text',label=paste0(' ','~R^2 == ',round(summ$r.squared,2)),
             parse=TRUE,x=Inf,y=-Inf,size=7*5/14, hjust=1, vjust=-2) +
    geom_text(aes(label=selectCS,vjust=selectCSv,color=selectCS),size=7*5/14) +
    labs(x='Crossover density (Predicted)', y='Crossover density (Observed)') +
    scale_color_manual(values = c('black','orange','red'))+
    scale_fill_manual(values = c('black','orange','red'))+
    theme(legend.position = 'none')
  
  #ggtitle(paste(as.character(f)[2],as.character(f)[1],as.character(f)[3],collapse = " "))+
  
  gRes <- ggplot(dfPlot,aes(x=chromSize,y=residual)) + geom_point() + geom_hline(yintercept=0)
  
  return(list(figFit=g,resid=gRes,fig=ggarrange(g,gRes,ncol=1,nrow=2),r2=summ$r.squared))
}

checkChroms <- function(df, p){
  ret <- list()
  ret$df <- data.frame(chrom=c("all",paste0("chr",1:19)),r2=0)
  
  ret$all   <- plotObsVsPred(df, res="COperMb", preds=p)  
  ret$df$r2[ret$df$chrom=="all"] <- ret$all$r2
  
  for (cs in 1:19){
    ret[[paste0("chr",cs)]]  <- plotObsVsPred(df[df$chr != paste0("chr",cs),], res="COperMb", preds=p)
    ret$df$r2[ret$df$chrom==paste0("chr",cs)] <- ret[[paste0("chr",cs)]]$r2
  }
  return (ret)
}

checkPredictors <- function(dfPlot, res="COperMb", preds=preds){
  
  dfRet <- data.frame(p=preds,r2=0)
  for (p in preds){
    print(p)
    f <- as.formula( paste(res, paste(p, collapse = " + "), sep = " ~ "))
    
    fit <- eval(bquote(lm(.(f), data = dfPlot )) )
    summ <- summary(fit)
    dfRet$R2[dfRet$p == p] <- summ$r.squared
  }
  
  return(dfRet)
}

predictCOperMb <- function(repSpeed){
  
  ## Get sizes
  x         <- fread('Gapped.mm10.bed', col.names = c("chr","start","end"))
  mm10Chrom <- x[ , .(chromSize = sum(end - start)), by = .(chr) ]
  mm10Chrom <- mm10Chrom[mm10Chrom$chr %in% paste0('chr',1:19),]
  
  ## Get other data
  CO        <- read.table ('B6CAST_Crossoverss.bed', col.names = c("chr","start","end"), stringsAsFactors = FALSE)
  DSB       <- read.table ('B6fXCASTm_hotspots.bedgraph',col.names = c("chr","start","end","str"), stringsAsFactors = FALSE)
  ori       <- read.table('hiconf_origins.mm10.bedgraph',col.names = c("chr","start","end","str"), stringsAsFactors = FALSE)
  gc        <- read.table('cs.gc.bed', header = F, stringsAsFactors = FALSE); 
  names(gc) <- c('chr','gc')
  gc        <- gc[gc$chr %in% paste0('chr',1:19),]
  
  ## reshape the RT speed data
  rtSpeed <- as.data.frame(repSpeed %>% group_by(cs,pc) %>% dplyr::summarise(rt=round(mean(time))))
  
  head(rtSpeed)
  dfRTSpeed           <- rtSpeed[rtSpeed $pc == 10,c('cs','rt')]
  dfRTSpeed$chr       <- paste0('chr',dfRTSpeed$cs)
  names(dfRTSpeed)[2] <- 'cyc10'
  
  dfRTSpeed$cyc25 <- rtSpeed$rt[rtSpeed $pc == 25]
  dfRTSpeed$cyc50 <- rtSpeed$rt[rtSpeed $pc == 50]
  dfRTSpeed$cyc75 <- rtSpeed$rt[rtSpeed $pc == 75]
  dfRTSpeed$cyc90 <- rtSpeed$rt[rtSpeed $pc == 90]
  
  ## Merge the different data into one df
  d1 <- data.frame(DSB %>% 
                     dplyr::filter(chr %in% paste0('chr',1:19)) %>% 
                     inner_join(mm10Chrom, by = 'chr') %>% 
                     group_by(chr) %>% 
                     dplyr::mutate(perMb=sum(str)/chromSize*1000000) %>% 
                     dplyr::summarise(DSBperMb=mean(perMb)))
  
  d2 <- data.frame(DSB %>% 
                     dplyr::filter(chr %in% paste0('chr',1:19)) %>% 
                     inner_join(mm10Chrom, by = 'chr') %>% 
                     group_by(chr) %>% 
                     dplyr::mutate(perMb=n()/chromSize*1000000) %>% 
                     dplyr::summarise(DSBHSperMb=mean(perMb)))
  
  d3 <- data.frame(ori %>% 
                     dplyr::filter(chr %in% paste0('chr',1:19)) %>% 
                     inner_join(mm10Chrom, by = 'chr') %>% 
                     group_by(chr) %>% 
                     dplyr::mutate(perMb=sum(str)/chromSize*1000000) %>% 
                     dplyr::summarise(OriperMb=mean(perMb)))
  
  d4 <- data.frame(ori %>% 
                     dplyr::filter(chr %in% paste0('chr',1:19)) %>% 
                     inner_join(mm10Chrom, by = 'chr') %>% 
                     group_by(chr) %>% 
                     dplyr::mutate(perMb=n()/chromSize*1000000) %>% 
                     dplyr::summarise(OriHSperMb=mean(perMb)))
  
  d5 <- data.frame(CO %>% 
                     dplyr::filter(chr %in% paste0('chr',1:19)) %>%
                     inner_join(mm10Chrom, by = 'chr') %>% 
                     group_by(chr) %>% 
                     dplyr::mutate(perMb=n()/chromSize*1000000) %>% 
                     dplyr::summarise(COperMb=mean(perMb)))
  
  all <- d1 %>% inner_join(d2, by = 'chr')
  all <- all %>% inner_join(d3, by = 'chr')
  all <- all %>% inner_join(d4, by = 'chr')
  all <- all %>% inner_join(d5, by = 'chr')
  all <- all %>% inner_join(gc, by = 'chr')
  all <- all %>% inner_join(dfRTSpeed, by = 'chr')
  all <- all %>% inner_join(mm10Chrom, by = 'chr')
  
  all$chr <- factor(all$chr,levels=paste0("chr",1:19))
  all$chromSize <- all$chromSize/1000000
  
  formula <- y ~ x
  
  #v1 <- checkChroms(all, p=c('poly(chromSize,2)'))
  
  ## For main figure; 
  dSize   <- plotObsVsPred(all, res="COperMb", preds=c('poly(chromSize,2)'), name = "Chr. size only")
  dRep    <- plotObsVsPred(all, res="COperMb", preds=c('cyc75'), name = "Replication speed only")
  dSzRep  <- plotObsVsPred(all, res="COperMb", preds=c('poly(chromSize,2) + cyc75'), name = "Size & Replication speed")
  dAll    <- plotObsVsPred(all, res="COperMb", preds=c('poly(chromSize,2) + cyc75 + DSBperMb'), name = "Size & Replication speed & DSBs")
  
  mainFigs <- list(sz=dSize$figFit,
                   rt=dRep$figFit,
                   x2=dSzRep$figFit,
                   all=dAll$figFit)
  
  ## For supp
  pI       <- checkPredictors(all,"COperMb",preds = c('chromSize','poly(chromSize,2)','gc','DSBperMb','DSBHSperMb','OriperMb','cyc10','cyc25','cyc50','cyc75','cyc90'))
  
  pI$nm[pI$p == 'poly(chromSize,2)'] <- 'Chr. size (poly)'
  pI$nm[pI$p == 'chromSize']         <- 'Chr. size (linear)'
  pI$nm[pI$p == 'gc']                <- 'GC (%)'
  pI$nm[pI$p == 'DSBperMb']          <- 'DSBs (per Mb)'
  pI$nm[pI$p == 'DSBHSperMb']        <- 'DSB hotspots (per Mb)'
  pI$nm[pI$p == 'OriperMb']          <- 'Origins (per Mb)'
  pI$nm[pI$p == 'cyc10']             <- 'Rep. Speed (10%)'
  pI$nm[pI$p == 'cyc25']             <- 'Rep. Speed (25%)'
  pI$nm[pI$p == 'cyc50']             <- 'Rep. Speed (50%)'
  pI$nm[pI$p == 'cyc75']             <- 'Rep. Speed (75%)'
  pI$nm[pI$p == 'cyc90']             <- 'Rep. Speed (90%)'
  
  pI$nm     <- factor(pI$nm,levels=pI$nm[order(pI$R2)])
  pI$barcol <- 'grey'
  
  gSuppI    <- ggplot(pI,aes(x=nm,y=R2,fill=barcol)) + 
    geom_bar(stat='identity',width = .8,lwd=.3,color='grey50') + 
    coord_flip() + 
    xlab('Chromosome property') + 
    ylab(bquote('Model '*R^2)) +
    theme(legend.position='none')
  #ylab(bquote('Model regression coefficient ('*R^2*')'))
  
  ## LMs with 2x predictors
  p2       <- checkPredictors(all,"COperMb",preds = c('poly(chromSize,2)', 
                                                      'poly(chromSize,2) + gc',
                                                      'poly(chromSize,2) + DSBperMb',
                                                      'poly(chromSize,2) + DSBHSperMb',
                                                      'poly(chromSize,2) + OriperMb',
                                                      'poly(chromSize,2) + cyc10',
                                                      'poly(chromSize,2) + cyc25',
                                                      'poly(chromSize,2) + cyc50',
                                                      'poly(chromSize,2) + cyc75',
                                                      'poly(chromSize,2) + cyc90'))
  
  p2$R2ref <- (p2$R2 - p2$R2[p2$p == 'poly(chromSize,2)'])/p2$R2[p2$p == 'poly(chromSize,2)']*100
  p2$p     <- factor(p2$p,levels=p2$p[order(p2$R2ref)])
  
  p2$nm = gsub(pattern = 'poly\\(chromSize,2\\)\\ \\+\\ ',replacement = '',p2$p)
  p2$nm[p2$nm == 'poly(chromSize,2)'] <- 'Chr. size (poly)'
  p2$nm[p2$nm == 'gc']                <- 'GC (%)'
  p2$nm[p2$nm == 'DSBperMb']          <- 'DSBs (per Mb)'
  p2$nm[p2$nm == 'DSBHSperMb']        <- 'DSB hotspots (per Mb)'
  p2$nm[p2$nm == 'OriperMb']          <- 'Origins (per Mb)'
  p2$nm[p2$nm == 'cyc10']             <- 'Rep. Speed (10%)'
  p2$nm[p2$nm == 'cyc25']             <- 'Rep. Speed (25%)'
  p2$nm[p2$nm == 'cyc50']             <- 'Rep. Speed (50%)'
  p2$nm[p2$nm == 'cyc75']             <- 'Rep. Speed (75%)'
  p2$nm[p2$nm == 'cyc90']             <- 'Rep. Speed (90%)'
  
  p2$nm     <- factor(p2$nm,levels=p2$nm[order(p2$R2ref)])
  p2$barcol <- 'grey'
  
  g2        <- ggplot(p2,aes(x=nm,y=R2,fill=barcol)) + 
    geom_bar(stat='identity',width = .8,lwd=.3,color='grey50') + 
    coord_flip() + 
    xlab('Chr. size (poly) + ') + 
    ylab(bquote('Model '*R^2)) + 
    theme(legend.position='none')
  #ylab(bquote('Model regression coefficient ('*R^2*')'))
  
  g2I       <- ggplot(p2) + geom_tile(aes(x=1, y=nm,fill=R2ref),color='grey40') + 
    geom_text(aes(x=1, y=nm, label=paste0(round(R2ref,0),'%')),size=7*5/14) + 
    scale_fill_gradient(low='white',high='red') + 
    ylab('Model improvement (%)') + 
    theme(axis.text = element_blank(), axis.ticks= element_blank(),axis.line.y=element_blank(),legend.position='none') + 
    xlab('')
  
  gSuppHoriz <- ggarrange(gSuppI,
                          g2,
                          g2I + theme(axis.text.y=element_blank()),
                          ncol=3,nrow=1,
                          widths=c(5,5,1),
                          align='h')
  
  gSuppMulti <- ggarrange(g2,
                          g2I + theme(axis.text.y=element_blank()),
                          ncol=2,nrow=1,
                          widths=c(5,1),
                          align='h')
  
  return(list(df=all, 
              dfSupp1=pI, 
              dfSupp2=p2, 
              main=mainFigs, 
              gSupp1=gSuppI, 
              gSupp2=gSuppMulti,
              gSuppHZ=gSuppHoriz,
              gSuppPart1=gSuppI, 
              gSuppPart2=g2, 
              gSuppPart3=g2I))
}

plotOne <- function(dfPred){
  g<-ggplot(dfPred, aes(y = prediction, x = COperMb)) +
    geom_smooth(method = "lm", se=T,alpha=0.25, lwd=.4, colour='grey') +
    geom_point(aes(fill = selectCS, color = selectCS), size=0.8) +
    facet_wrap(~type,ncol=4) +
    stat_regline_equation(label.x.npc='left',
                          label.y.npc='top',
                          vjust=1,
                          aes(label =  paste(..rr.label.., sep = "~~~~")),size=7*5/14,color='black') +
    labs(y='Crossover density (Predicted)', x='Crossover density (Observed)') +
    scale_color_manual('',values = c('orange','red','black'))+
    scale_fill_manual('',values = c('orange','red','black'))+
    theme(legend.position = c(1,0),
          legend.justification=c(1,0),
          legend.background=element_blank(),
          legend.key.size=unit(0.2,'cm'),
          strip.background=element_blank(),
          strip.text=element_text(size=7,face='bold',hjust=0),
          panel.border=element_rect(size=.2,fill=NA)) + 
    coord_cartesian(xlim=c(19.7,39.5),expand=FALSE,clip='off') + 
    scale_x_continuous(breaks=c(25,35)) + 
    scale_y_continuous(breaks=c(25,35)) 
  return(g)
}
#  OK make figures ... ######################################################################################
repSpeedData <- getReplicationSpeed(modelFile = 'rtSimModel.Rdata')

lmData  <- predictCOperMb(repSpeedData$df)

nMargin <- theme(plot.margin = unit(c(0,0,0,.1),'cm'))
#############################################################################################################

dfPredSz <- lmData$main$sz$data %>% select(chr,COperMb,prediction,selectCS) %>% mutate(type='Chromosome Size',selectCS=ifelse(selectCS == '','others',selectCS))
dfRTSz <- lmData$main$rt$data   %>% select(chr,COperMb,prediction,selectCS) %>% mutate(type='Replication Speed',selectCS=ifelse(selectCS == '','others',selectCS))
dfX2Sz <- lmData$main$x2$data   %>% select(chr,COperMb,prediction,selectCS) %>% mutate(type='Size + Rep. Speed',selectCS=ifelse(selectCS == '','others',selectCS))
dfAllSz <- lmData$main$all$data %>% select(chr,COperMb,prediction,selectCS) %>% mutate(type='Size + Rep. Speed + DSBs',selectCS=ifelse(selectCS == '','others',selectCS))

df <- rbind(dfPredSz,dfRTSz,dfX2Sz,dfAllSz)

noLeg <- theme(legend.position='none',axis.text.y=element_blank(),axis.ticks.y=element_blank())

gA <- plotOne(dfPredSz) + noLeg + nMargin + coord_cartesian(xlim=c(19.7,39.5))
gC <- plotOne(dfRTSz) + noLeg + nMargin + ylab('')+ coord_cartesian(xlim=c(19.7,39.5))
gF <- plotOne(dfX2Sz) + noLeg + nMargin + ylab('')+ coord_cartesian(xlim=c(19.7,39.5))
gG <- plotOne(dfAllSz) + nMargin + ylab('')+ coord_cartesian(xlim=c(19.7,39.5))

gPreds <- ggarrange(gA,gC,gF,gG,ncol=4,
                    labels=c('A','C','F','G'),
                    font.label = list(size=8,fontface='bold'))

#gPreds <- annotate_figure(gPreds, bottom = text_grob(label = 'Crossover density (observed)',size=7))

dfRepSpeed <- repSpeedData$df %>% group_by(cs,pc,.drop = FALSE) %>% 
  summarise(median=median(time),pc5=quantile(time,.05),pc95=quantile(time,.95)) %>% 
  dplyr::filter(pc %in% c(25,50,75)) %>%
  mutate(pc=paste0(pc,'%'),selectCS=ifelse(cs==11,'chr11',ifelse(cs==13,'chr13','other')))

gB <- ggplot(dfRepSpeed,aes(color=selectCS)) + 
  geom_segment(aes(x=pc5,xend=pc95,y=cs,yend=cs),lwd=.2) + 
  geom_segment(aes(x=median,xend=median,y=cs+0.6,yend=cs-0.6),lwd=.5) + 
  facet_wrap(~pc,ncol=1,switch = TRUE) + xlab('Time to replicate (RT-Sim cycles)') + 
  ylab('Chromosome') + 
  scale_y_continuous(breaks=c(1,10,19)) +
  scale_color_manual('',values = c('orange','red','black')) + 
  theme(legend.position = c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank(),
        legend.key.size=unit(0.2,'cm'),
        strip.background=element_rect(fill='grey90'),
        strip.text=element_text(size=7,face='bold',hjust=0.5)) 

nDEMargin <- theme(plot.margin = unit(c(0.4,0,0,0),'cm'))

gDE <- ggarrange(lmData$gSuppPart1 + nDEMargin + scale_fill_manual(values=c('grey70')) + scale_y_continuous(breaks=c(0,0.25,0.5,.75,1),labels=c(0,'','','',1)) + coord_flip(xlim=c(0.5,11.5),ylim=c(0,1.01),expand=FALSE),
                 lmData$gSuppPart2 + nDEMargin + scale_fill_manual(values=c('grey70')) + scale_y_continuous(breaks=c(0,0.25,0.5,.75,1),labels=c(0,'','','',1)) + coord_flip(xlim=c(0.5,10.5),ylim=c(0,1.01),expand=FALSE),
                 lmData$gSuppPart3 + nDEMargin + theme(axis.text.y=element_blank()),
                 ncol=3,nrow=1,
                 widths=c(5,5,1),
                 align='h',
                 labels = c('D','E',''),
                 vjust = 1,hjust=0,
                 font.label = list(size=8,face='bold',vjust=1)) 

#gCD <- annotate_figure(gCD, bottom = text_grob(label = bquote('Model correlation coefficient ('*R^2*')'),
#                                               size=7))

gBDE <- ggarrange(gB,gDE,
                  ncol=2,
                  widths=c(1,2),
                  labels = c('B',''),
                  vjust = 1,hjust=0,
                  font.label = list(size=8,face='bold',vjust=1))

gFig <- ggarrange(gPreds,
                   gBDE,
                   nrow=2,
                   heights=c(5,6),
                   labels = c('',''),
                   vjust = 1,hjust=0,
                   font.label = list(size=8,face='bold',vjust=1))

ggsave('Pratto_et_al_Figure6.png',gFig,height=3.5,width=6.5,dpi = 400)
ggsave('Pratto_et_al_Figure6.pdf',gFig,height=3.5,width=6.5)
