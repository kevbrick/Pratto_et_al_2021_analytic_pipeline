imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

load('hiconf_origins.mm10.Rdata')

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(data.table)

theme7point()

calc_auc <- function(X, Y){
  # inputs already sorted, best scores first 
  dY <- c(diff(Y), 0)
  dX <- c(diff(X), 0)
  sum(Y * dX) + sum(dY * dX)/2
}

plotABCF <- function(x){
  z <- read.table('closest.tab')
  names(z) <- c('pos','ori','type')

  z$pos <- as.numeric(as.character(z$pos))

  zT <- as.data.frame(table(z))

  zT$pos  <- as.numeric(as.character(zT$pos))
  zT$Freq <- as.numeric(as.character(zT$Freq))

  zB <-  by(zT$Freq,zT[,c('ori','type')],FUN=cumsum)

  zT$cumsum <- c(zB[[1]],zB[[2]],zB[[3]],zB[[4]],zB[[5]],zB[[6]])/max(zB[[1]])*100

  zP <- zT[zT$type == 'tss',]

  g2 <- ggplot(zP,aes(x=pos,y=cumsum,color=ori)) +
    geom_vline(xintercept=1000,lwd=.2,lty='dashed',color='grey50') +
    geom_hline(yintercept=66,  lwd=.2,lty='dashed',color='grey50') +
    geom_line() +
    scale_x_log10(breaks=c(1,1000,1000000),
                  labels=c('1 bp','1 Kbp', '1 Mbp')) +
    annotation_logticks(sides='b',
                        long=unit(.2,'cm'),
                        mid=unit(.1,'cm'),
                        short=unit(.1,'cm'),
                        size=.2) +
    xlab('Distance to closest TSS') +
    ylab('Cumulative Percentage (%)') +
    theme(legend.position = c(0,1),
          legend.justification=c(0,1),
          legend.title=element_blank(),
          legend.background=element_blank(),
          legend.key.size = unit(0.2,'cm')) +
    coord_cartesian(ylim=c(-5,100)) +
    annotate(geom='text',label='66 %',x=400,y=75,size=8*5/14) +
    scale_color_manual(values=c('forestgreen','grey50'))

  x <- read.table('oriVori.tab')

  g1 <- ggplot(x,aes(x=V1,color=V2,fill=V2)) + geom_density(alpha=.2,adjust=2) +
    scale_x_log10(breaks=c(1,1000,1000000),
                  labels=c('1 bp','1 Kbp', '1 Mbp')) +
    annotation_logticks(sides='b',
                        long=unit(.2,'cm'),
                        mid=unit(.1,'cm'),
                        short=unit(.1,'cm'),
                        size=.2) +
    xlab('Inter-origin distance') +
    ylab('Density') +
    theme(legend.position = c(0,1),
          legend.justification=c(0,1),
          legend.title=element_blank(),
          legend.background=element_blank()) +
    scale_color_manual(values=c('forestgreen','grey50')) +
    scale_fill_manual(values=c('forestgreen','grey50'))
  #coord_cartesian(ylim=c(-5,100))

  g4        <- read.table('g4_vs_origins.tab',header=FALSE)
  names(g4) <- c('ori','g4','d','strand')
  g4$d2c    <- floor(g4$d/50)*50

  dP        <- as.data.frame(table(g4[,c('d2c')]))
  dP$pos    <- as.numeric(as.character(dP$Var1))
  dP$N      <- as.numeric(as.character(dP$Freq))
  dP$FC     <- dP$N / mean(dP$N[abs(dP$pos) > 4000])

  g4plot <- ggplot(dP,aes(x=pos/1000,y=FC)) +
    geom_vline(xintercept=0,lwd=.2,lty='dashed',color='grey50') +
    geom_hline(yintercept=1,lwd=.2,color='grey50') +
    geom_smooth(span=.1,lwd=.4) +
    coord_cartesian(xlim=c(-5,5),expand=FALSE) +
    geom_point(size=.2) +
    xlab('Distance to origin center (Kb)') +
    ylab('G4 density\n(fold-change / flanks)') +
    scale_x_continuous(breaks=c(-4,-2,0,2,4))

  ################################################################### BENDABILITY
  
  bData        <- fread('mm10Origins.bendability.txt',header=FALSE,drop=4)
  #bData        <- fread('mb.txt',header=FALSE,drop=4)
  names(bData) <- c('seq','pos','bend')
  bData$name   <- 'Origins'
  
  ## Add randoms
  bR        <- bData
  bR$pos    <- sample(bR$pos,replace=FALSE,size=length(bR$pos))
  bR$name   <- 'Random'
  
  df <- rbind(bData,bR)
  
  iD <- df %>%
    group_by(pos,name) %>%
    dplyr::summarise(mean=mean(bend), sd=sd(bend))
  
  gMN <- ggplot(iD,aes(x=(10000-pos*3)/1000,y=mean,color=name,fill=name)) + 
    geom_point(size=.2,alpha=.8,show.legend = FALSE) + 
    geom_smooth(span=.1,lwd=.3) + 
    xlab('Distance to origin center (Kb)') + 
    ylab('DNA bendability (Mean)') +
    scale_color_manual(values=c('forestgreen','grey70')) + 
    scale_fill_manual(values=c('palegreen3','grey80')) + 
    theme(legend.title=element_blank(),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.direction = 'vertical',
          legend.key.size=unit(0.2,'cm')) + 
    coord_cartesian(xlim=c(-5,5),expand=FALSE) +
    scale_x_continuous(breaks=c(-4,-2,0,2,4))
  
  gSD <- ggplot(iD,aes(x=(10000-pos*3)/1000,y=sd  ,color=name,fill=name)) + 
    geom_point(size=.2,alpha=.8) + 
    geom_smooth(span=.1,lwd=.3) + 
    xlab('Distance to origin center (Kb)') + 
    ylab('DNA bendability (S.D.)') +
    scale_color_manual(values=c('forestgreen','grey70')) + 
    scale_fill_manual(values=c('palegreen3','grey80')) + 
    theme(legend.title=element_blank(),
          legend.position='none') + 
    coord_cartesian(xlim=c(-5,5),expand=FALSE) +
    scale_x_continuous(breaks=c(-4,-2,0,2,4))
  
  gBendability <- ggarrange(gMN + xlab(''),gSD,ncol=1,nrow=2,align='v',heights=c(2,1))
  ####################################################################
  
  gPanelAB <- ggarrange(g1,g2,
                   ncol=1,nrow=2,
                   labels=c('A','B'),
                   align='v',
                   vjust=1,
                   font.label = list(size = 8,
                                     face='bold'))

  gPanelCF <- ggarrange(gMN + xlab(''),gSD,g4plot,
                        ncol=1,nrow=3,
                        labels=c('C','','F'),
                        align='v',
                        vjust=1,
                        heights=c(2,1,2),
                        font.label = list(size = 8,
                                          face='bold'))
  
  gPanelABCF <- ggarrange(gPanelAB,gPanelCF,
                        ncol=2,nrow=1,
                        align='h',
                        vjust=1,
                        font.label = list(size = 8,
                                          face='bold'))
  return(list(gA=g1, gB=g2, gC = gBendability, gF = g4plot, gAll =gPanelABCF));
}

makeBarPlots <- function(){
  
  ## Barplot Data:
  dfATACCGI <- read.table('figure1eData.tab',header=TRUE) 
  
  dfPlotOri <- dfATACCGI %>% mutate(pc=N/tot*100,
                                    lbl=paste0(" ",round(N/tot*100,0),'%'))
  
  names(dfPlotOri)[3] <- 'ATAC\nSeq'
  
  dfPlotOri$name <- factor(dfPlotOri$name, 
                           levels=rev(c('ATAC-CGI','ATACnoCGI','CGInoATAC','Neither')))
  
  gOriBars <- ggplot(dfPlotOri,aes(x=pc,y=name)) + 
    geom_bar(stat='identity',width = .6, fill='grey50',color='black',lwd=.3) + 
    xlab('Origins (%)') + 
    coord_cartesian(xlim=c(0,90),ylim=c(0.5,4.5),expand=FALSE) +
    theme(legend.position='none',
          axis.text.y=element_blank(),
          plot.margin=unit(c(0,0,0,0),'cm')) + 
    geom_text(aes(label=lbl),size=7*5/14,hjust=-0.5,color='black')+
    ylab('')
  
  dfKeyOri <- reshape2:::melt.data.frame(dfPlotOri,id.vars = c('name','pc'),
                                         measure.vars = c('ATAC\nSeq','CGI'))
  
  gKeyOri <- ggplot(dfKeyOri,aes(x=variable,y=name,fill=value)) + 
    geom_tile(color=NA,fill=NA,lwd=.2) + 
    scale_fill_manual(values=c('white','grey90')) + 
    geom_text(size=7*5/14,aes(label=ifelse(value,'+','-'))) + 
    theme(legend.position='none',
          axis.text.y=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          plot.margin=unit(c(0,0,0,0),'cm')) + 
    coord_cartesian(ylim=c(0.5,4.5),expand=FALSE) +
    xlab('') + ylab('')
  
  gOri <- ggarrange(gKeyOri,
                    gOriBars,
                    ncol=2,nrow=1,
                    widths=c(1,2),
                    align='h')
  
  ## Figure for sites at origins -----------------------------------------------
  dfATACCGI <- fread('figure1fData.tab',header=TRUE) 
  dfPlotSites <- dfATACCGI %>% mutate(pc=N/tot*100,
                                      lbl=paste0("  ",round(N/tot*100,0),'%'),
                                      type=ifelse(type=='ATAC-CGI','ATAC & CGI',type))
  
  dfPlotSites$type <- factor(dfPlotSites$type, 
                             levels=rev(c('ATAC & CGI','CGI no ATAC','ATAC no CGI')))
  
  gSitesAtOriInit <- ggplot(dfPlotSites,aes(x=pc,y=type)) + 
    geom_bar(stat='identity',width = .5, fill='grey50',color='black',lwd=.3) + 
    xlab('At origins (%)') + 
    coord_cartesian(xlim=c(0,100),ylim=c(0.5,3.5),clip='off',expand=FALSE) +
    theme(legend.position='none',
          plot.margin=unit(c(0,0.3,0,0),'cm')) + 
    geom_text(aes(label=lbl),size=7*5/14,hjust=-0.5,color='black')+
    ylab('')
  
  gSitesAtOri <- ggarrange(ggplot() + theme_void(),
                           gSitesAtOriInit,
                           ncol=2,nrow=1,
                           widths=c(1,7),
                           align='h')
  
  gX2 <- ggarrange(gOriBars,gSitesAtOri,
                   ncol=1,nrow=2,
                   heights=c(5,5),
                   labels = c('G','H'),
                   vjust = 1,hjust=0,
                   font.label = list(size=8,face='bold',vjust=1))
  
  return(list(gOrisWith=gOri,
              gLociWith=gSitesAtOriInit,
              gX2=gX2))
}

plotCpGpredictiveAbility <- function(figLbl1='A',figLbl2='B'){
  dfROCCpG <- read.table('CpG_overlapStats.txt',header=TRUE)
  
  dfROCCpG <- dfROCCpG[order(dfROCCpG$pcGC,decreasing = TRUE),]
  
  #oops ... this is PPV ## 
  dfROCCpG$fp           <- dfROCCpG$totGC-dfROCCpG$atOri
  dfROCCpG$sens         <- dfROCCpG$atOri/dfROCCpG$totOri
  dfROCCpG$spec         <- (max(dfROCCpG$fp)-dfROCCpG$fp)/max(dfROCCpG$fp)
  dfROCCpG$oneMinusSpec <- 1 - dfROCCpG$spec
  
  AUC <- round(calc_auc(dfROCCpG$oneMinusSpec,dfROCCpG$sens),2)
  
  dfROCTF <- read.table('ROCdata.tab',header=TRUE)
  
  dfROC <- rbind(dfROCTF %>% mutate(type='TF:TFBS') %>% select(type,name,spec,sens),
                dfROCCpG %>% mutate(type="Ori:CpG",name="OriCpG",spec,sens)%>% select(type,name,spec,sens))
  
  dfROC$type <- factor(dfROC$type,levels=c('Ori:CpG','TF:TFBS'))
  gROC <- ggplot() + 
    geom_line(data=dfROC ,aes(x=1-spec,y=sens,color=type,group=name),lwd=.2) + 
    geom_line(data=data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color='grey60',lwd=.2) + 
    geom_point(data=dfROC[dfROC$type=='Ori:CpG',],aes(x=1-spec,y=sens,color=type),size=.5) + 
    geom_line(data=dfROC[dfROC$type=='Ori:CpG',] ,aes(x=1-spec,y=sens,color=type,group=name),lwd=.4) + 
    coord_cartesian(xlim=c(-0.02,1.025),ylim=c(-0.02,1.025),expand=FALSE) + 
    xlab('1 - Specificity') + 
    ylab('Sensitivity') + 
    scale_color_manual('',values=c('forestgreen','grey70')) + 
    annotate(geom='text',label=paste0("AUC = ",AUC),x=.5,y=1,hjust=0.5,size=7*5/14,color='forestgreen') + 
    theme(legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background = element_blank(),
          legend.key.size=unit(0.3,'cm'))
  
  ms <- reshape2:::melt.data.frame(dfROCCpG,id.vars = c('pcGC','spec'), measure.vars=c('atOri','fp'))
  ms$type <- "False positive"
  ms$type[ms$variable == 'atOri'] <- "At origin"
  
  sd <- dfROCCpG[dfROCCpG$pcGC == 0.07,]
  lbl <- paste0(round(sd$sens*100,0),'% of origins\n(N = ',sd$atOri,'; FP = ',sd$fp,')')
  
  gOri <- ggplot(ms,aes(x=pcGC*100,y=value)) + 
    geom_vline(xintercept=7,lwd=.2,color='magenta',lty='dashed')+
    geom_hline(yintercept=ms$value[ms$pcGC == 0.07 & ms$variable == 'atOri'],lwd=.2,color='magenta',lty='dashed')+
    geom_point(size=.4,aes(color=type)) + geom_line(aes(color=type),lwd=.2) + 
    scale_y_log10(labels=fancy_scientific) +
    scale_x_reverse() +
    xlab('GC content threshold (%)') + 
    ylab('CpG peaks (#)') +
    annotation_logticks(sides='l',
                        short=unit(0.1,'cm'),
                        mid=unit(0.1,'cm'),
                        long=unit(0.1,'cm'),
                        size=.2) + 
    scale_color_manual(values=c('forestgreen','grey60')) +
    theme(legend.title=element_blank(),
          legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background=element_blank()) + 
    annotate(geom='text',x=15,y=50000,
             hjust=0,
             label=lbl,color='magenta',size=7*5/14)
  
  gX2 <- ggarrange(gOri,gROC,
                   ncol=1,nrow=2,
                   labels=c(figLbl1,figLbl2),
                   align='v',
                   vjust=1,
                   heights=c(1,1),
                   font.label = list(size = 8,
                                     face='bold'))
  return(list(gRaw = gOri, 
              gROC = gROC,
              gX2 = gX2))
}

load('hiconf_origins.mm10.Rdata')
diNT <- plotDiNT_v_OrisMM10(originsDF)

##########
gABCFdets <- plotABCF()
gIJ       <- plotCpGpredictiveAbility('G','H')

marg <- theme(plot.margin=unit(c(0,0,0,0.4),"cm"))
gDi  <- ggarrange(diNT$gCenter + theme(plot.margin=unit(c(0.4,0,0,0.4),"cm")),
                  diNT$gCC + theme(plot.margin=unit(c(0,0,0,0.4),"cm")),
                  ncol=1,nrow=2,
                  labels=c('D','E'),
                  vjust=1,
                  font.label = list(size = 8,
                                    face='bold'))

gBoth  <- ggarrange(gABCFdets$gA,gDi,
                  ncol=1,nrow=2,
                  heights=c(7,5))

lstBars <- makeBarPlots()

################################################################################
gAB <- ggarrange(gABCFdets$gA,
                  gABCFdets$gB,
                  ncol=2,nrow=1,
                  labels=c('A','B'),
                  vjust=1,
                  font.label = list(size = 8,
                                    face='bold'))

gABDE <- ggarrange(gAB,gDi,
                   ncol=1,nrow=2,
                   heights=c(2,3))

gCGH <- ggarrange(gABCFdets$gC,
                  lstBars$gOrisWith,
                  lstBars$gLociWith,
                  ncol=1,nrow=3,
                  heights=c(3.5,2,2),
                  labels=c('C','G','H'),
                  vjust=1,
                  font.label = list(size = 8, face='bold'))

gABCDEGH <- ggarrange(gABDE,gCGH,
                      nrow=1,ncol=2,
                      widths=c(2,1))

gFIJ <- ggarrange(gABCFdets$gF,
                  gIJ$gRaw,
                  gIJ$gROC,
                  ncol=3,nrow=1,
                  labels=c('F','I','J'),
                  vjust=1,
                  font.label = list(size = 8,
                                    face='bold'))

gBoth <- ggarrange(gABCDEGH,gFIJ,
                   ncol=1,nrow=2,heights=c(5,2))
nScale <- 1
ggsave(plot = gBoth,
       filename = 'Pratto_et_al_SuppFig_origin_properties_and_CpGs.png',
       width=6.5*nScale, height=8*nScale)

ggsave(plot = gBoth,
       filename = 'Pratto_et_al_SuppFig_origin_properties_and_CpGs.pdf',
       width=6.5*nScale, height=8*nScale)


