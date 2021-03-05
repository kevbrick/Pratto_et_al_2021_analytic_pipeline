source('accessoryFiles/scripts/R/repliSim_loadModules.R')
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

library(data.table)

theme7point()

colorBlind8   <- c("SpgB" = "#009E73", 
                   "SpgU"="#CC79A7", 
                   "SpgI"="#D55E00", 
                   "MeiS" = "#444444",
                   "Unclassified" = "#cccccc", 
                   "Other1" = "#E69F00", 
                   "Other2" = "#56B4E9", 
                   "Other3" = "#0072B2",
                   "Other4" = "#F0E442")

plotSorting <- function(colz    = colorBlind8,
                        pointSz = 0.1,
                        pointA  = 0.2){

  theme7point()

  sd1 <- getSort1Data()

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

  dRep$stage <- "Unclassified"
  dRep$stage[dRep$DMRT1 > 3000 & dRep$DMRT1 < 11000 & dRep$STRA8 < 2000  & dRep$STRA8 > 550 ] <- 'SpgU'
  dRep$stage[dRep$DMRT1 > 7000 & dRep$DMRT1 < 30000 & dRep$STRA8 < 30000  & dRep$STRA8 > 2500] <- 'SpgB'
  dRep$stage[dRep$DMRT1 < 3000 & dRep$DMRT1 > 200 & dRep$STRA8 < 30000  & dRep$STRA8 > 6000] <- 'MeiS'
  dRep$stage[dRep$DMRT1 > 11000 & dRep$DMRT1 < 30000 & dRep$STRA8 < 2000  & dRep$STRA8 > 550] <- 'SpgI'

  d4C$stage <- "Unclassified"
  d4C$stage[d4C$DMRT1 > 3000 & d4C$DMRT1 < 11000 & d4C$STRA8 < 2000  & d4C$STRA8 > 550 ] <- 'SpgU'
  d4C$stage[d4C$DMRT1 > 7000 & d4C$DMRT1 < 30000 & d4C$STRA8 < 30000  & d4C$STRA8 > 2500] <- 'SpgB'
  d4C$stage[d4C$DMRT1 < 3000 & d4C$DMRT1 > 200 & d4C$STRA8 < 30000  & d4C$STRA8 > 6000] <- 'MeiS'
  d4C$stage[d4C$DMRT1 > 11000 & d4C$DMRT1 < 30000 & d4C$STRA8 < 2000  & d4C$STRA8 > 550] <- 'SpgI'

  d2C$stage <- "Unclassified"
  d2C$stage[d2C$DMRT1 > 3000 & d2C$DMRT1 < 11000 & d2C$STRA8 < 2000  & d2C$STRA8 > 550 ] <- 'SpgU'
  d2C$stage[d2C$DMRT1 > 7000 & d2C$DMRT1 < 30000 & d2C$STRA8 < 30000  & d2C$STRA8 > 2500] <- 'SpgB'
  d2C$stage[d2C$DMRT1 < 3000 & d2C$DMRT1 > 200 & d2C$STRA8 < 30000  & d2C$STRA8 > 6000] <- 'MeiS'
  d2C$stage[d2C$DMRT1 > 11000 & d2C$DMRT1 < 30000 & d2C$STRA8 < 2000  & d2C$STRA8 > 550] <- 'SpgI'

  dRep$stage <- factor(dRep$stage, levels=c('MeiS','SpgU','SpgB','SpgI','Unclassified'))

  fortxt <- data.frame(stage = c('MeiS','SpgU','SpgB','SpgI'),
                       STRA8 = c(15000, 1000, 15000,600),
                       DMRT1 = c(2000, 3500, 20000,20000))

  gSort <-  ggplot(dRep,aes(x=STRA8,y=DMRT1,
                            fill=as.factor(stage),
                            color=as.factor(stage))) +
    geom_density2d(lwd=.05,alpha=.8) +
    geom_point(shape=20,size=pointSz,alpha=pointA) +
    scale_x_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) +
    scale_y_log10(breaks=c(10,100,1000,10000,100000), labels=fancy_scientific) +
    annotate(geom='rect',
             xmax=30000,xmin=6000,
             ymax=3000,ymin=200,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'MeiS'])+
    annotate(geom='rect',
             xmax=2000,xmin=550,
             ymax=11000,ymin=3000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'SpgU'])+
    annotate(geom='rect',
             xmax=30000,xmin=2500,
             ymax=30000,ymin=7000,
             lwd=.3, lty='dashed',
             fill=alpha('white',0),color=colz[names(colz) == 'SpgB'])+
    annotate(geom='rect',
             xmax=2000,xmin=550,
             ymax=30000,ymin=11000,
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
    coord_cartesian(ylim=c(20,50000),xlim=c(100,50000))+
    annotate(geom='text',
             x=30000,y=200,
             label='MeiS',
             hjust=1,vjust=1.5,
             size=8*5/14,
             color=colz[names(colz) == 'MeiS']) +
    annotate(geom='text',
              x=550,y=3000,
              label='SpgU',
              hjust=0,vjust=1.5,
              size=8*5/14,
              color=colz[names(colz) == 'SpgU'])+
    annotate(geom='text',x=30000,y=30000,
             label='SpgB',
             hjust=1,vjust=-.5,
             size=8*5/14,
             color=colz[names(colz) == 'SpgB']) +
  annotate(geom='text',
              x=550,y=30000,
              label='SpgI',
              hjust=0,vjust=-0.5,
              size=8*5/14,
              color=colz[names(colz) == 'SpgI'])

    gDens <- ggplot(dRep[!dRep$stage=="Unclassified",], aes(SYCP3, color=stage))+ 
      geom_density(size=0.5)+
      scale_color_manual(values=colz) + 
      xlab('SYCP3 fluorescence')+
    theme(legend.title=element_blank(),
          legend.position='none',
          legend.key.height = unit(.4,'cm'),
          legend.justification=c(1,1),
          legend.direction = 'vertical',
          legend.key.size = unit(.5,'cm'),
          legend.background=element_blank())+coord_cartesian(xlim=c(74.5,900)) +
      ylab('Relative density') +
      scale_y_continuous(breaks=c(0,0.005,0.010,0.015),labels=c(0,0.05,0.10,0.15)) +
      annotate(geom='text',
               x=750,y=Inf,
               label='MeiS',
               hjust=0,vjust=1,
               size=8*5/14,
               color=colz[names(colz) == 'MeiS']) +
      annotate(geom='text',
               x=750,y=Inf,
               label='SpgU',
               hjust=0,vjust=3,
               size=8*5/14,
               color=colz[names(colz) == 'SpgU'])+
      annotate(geom='text',
               x=750,y=Inf,
               label='SpgB',
               hjust=0,vjust=5,
               size=8*5/14,
               color=colz[names(colz) == 'SpgB'])+
      annotate(geom='text',
               x=750,y=Inf,
               label='SpgI',
               hjust=0,vjust=7,
               size=8*5/14,
               color=colz[names(colz) == 'SpgI'])

  return(list(gA = gSort,
              gB = gDens,
              gC = gDNAContent))
}

plotRT <- function(x,RT,
                   yLabel='',
                   colUp='#7fbf7b',
                   colDown='#af8dc3',
                   smoothSpan=100,
                   figSmooth=10,
                   csNorm=FALSE){
  #############################################
  ## Build DF for plot (speeds it up ... )
  pDF <- data.frame(cs=x$cs,
                    pos=(x$from+x$to)/2/1000000,
                    initRT=x[[RT]])

  if (csNorm){
    pDF$initRT <- pDF$iniTRT/mean(pDF$initRT)
  }

  plotDF <- data.frame(pos=rollapply(data = pDF$pos,
                                     smoothSpan,
                                     function(x){mean(x)}),
                       RT=rollapply(data = pDF$initRT,
                                    smoothSpan,
                                    function(x){mean(x)}),
                       cs = x$cs[1])


  #############################################
  ## Build RT up and RT down segments
  plotDF$RTup                  <- plotDF$RT
  plotDF$RTup[plotDF$RT < 0]   <- 0

  plotDF$RTdown                <- plotDF$RT
  plotDF$RTdown[plotDF$RT > 0] <- 0

  plotDF <- plotDF[seq(1,dim(plotDF)[1],by = figSmooth),]

  xMin <- 0
  xMax <- max(x$from)/1000000

  yRange <- max(abs(plotDF$RT))
  yMin   <- -1*yRange
  yMax   <- yRange
  yHi    <- roundToOneDigit(yMax*0.75,'down')
  yLo    <- roundToOneDigit(yMin*0.75,'down')

  if (yLabel == ''){
    # The double atop here is to minimize spacing between the two lines of the axis label
    yAxisLabel <- bquote(atop("",
                              atop("RT",
                              log[2]*"("*cov[n]*"/median)")))
  }else{
    # The double atop here is to minimize spacing between the two lines of the axis label
    yAxisLabel <- bquote(atop("",
                              atop("RT; "*.(yLabel),
                                  log[2]*"("*cov[n]*"/median)")))
  }

  #geom_line(aes(x=pos,y=RT),color='black',alpha=1,lwd=.3) +
  g <- ggplot(plotDF) +
    geom_hline(yintercept=0,color='grey10',alpha=.3,lwd=.2) +
    geom_area(aes(x=pos,y=RTup),fill=colUp,alpha=.8,color=darkenColor(colUp,1.5),lwd=.3) +
    geom_area(aes(x=pos,y=RTdown),fill=colDown,alpha=.8,color=darkenColor(colDown,1.5),lwd=.3) +
    geom_hline(yintercept=0,color='grey20',alpha=.6,lwd=.1) +
    ylab(yAxisLabel) + xlab('Position (Mb)') +
    scale_y_continuous(breaks=c(yLo,0,yHi)) +
    coord_cartesian(expand = FALSE,
                    xlim=c(xMin,xMax),
                    ylim=c(yMin,yMax)) + 
    theme(axis.title.y  = element_text(size = 10,margin = margin(t = -2, r = 0, b =0, l = -2, unit = "mm")),
          axis.title.x  = element_text(size = 10,margin = margin(t = -2, r = 0, b =0, l = 0, unit = "mm")))


  return(g)
}

plotCCHeatMap <- function(mDF){

  g <- ggCorMat(cor(mDF)^2,flipIt = TRUE,
                keepLeadingZeros = FALSE,
                decimalPlaces = 1,
                yOnRight = TRUE,
                noDiagonal = TRUE,
                tileFontScale = .8,
                xTilt=90) +
    scale_fill_gradient(low='white',high='red',na.value = NA) +
    annotate(geom='text',label='Spearman R',x=7,y=9,angle=45,size=7*5/14)

  return(g)
}

csNum <- 19
cs2use <- paste0('chr',csNum)

rr <- fread('rep_v_rec_MM10.tab',header=TRUE)

rNames <- names(rr);
n <- c(which(rNames == 'expRT_2C_ALL_S2_2C_NA_NA_mm10'),
       which(rNames == 'expRT_MeiS_ALL_S3_2to4C_STRA8_DMRT1_mm10'),
       which(rNames == 'expRT_MeiS_ALL_S4_2to4C_STRA8_DMRT1_mm10'),
       which(rNames == 'expRT_MeiS_ALL_S5_2to4C_STRA8_DMRT1_mm10'),
       which(rNames == 'expRT_SPGB_ALL_S2_2to4C_STRA8DMRT1hi_NA_mm10'),
       which(rNames == 'expRT_SPGB_ALL_S4_2to4C_STRA8DMRT1_NA_mm10'),
       which(rNames == 'expRT_SPGI_ALL_S1_2to4C_DMRT1hi_STRA8_mm10'),
       which(rNames == 'expRT_SPGI_ALL_S4_2to4C_DMRT1hi_STRA8_mm10'),
       which(rNames == 'expRT_SPGU_ALL_S2_2to4C_DMRT1_STRA8_mm10'),
       which(rNames == 'expRT_SPGU_ALL_S4_2to4C_DMRT1_STRA8_mm10'),
       which(rNames == 'expRT_ESC_ALL_S1_2to4C_NA_NA_mm10'),
       which(rNames == 'expRT_E14_Dey_R1_NA_NA_NA_mm10'),
       which(rNames == 'expRT_SSC_Yehuda_R1_NA_NA_NA_mm10'),
       which(rNames == 'expRT_PGC_Yehuda_R1_NA_NA_NA_mm10'),
       which(rNames == 'expRT_CD8_Yehuda_R1_NA_NA_NA_mm10'),
       which(rNames == 'expRT_J185a_Myoblast_Int61896107_mm10_mm10'),
       which(rNames == 'expRT_L1210_Lymphoblastoid_Ext73945012_mm10_mm10'),
       which(rNames == 'expRT_BCell_Tubbs_R1_NA_NA_NA_mm10'))

namez <- c('2C',
           'MeiS(a)',
           'MeiS(b)',
           'MeiS(c)',
           'SpgB(a)',
           'SpgB(b)',
           'SpgI(a)',
           'SpgI(b)',
           'SpgU(a)',
           'SpgU(b)',
           'ESC(a)',
           'ESC(b)',
           'SSC',
           'PGC',
           'CD8',
           'Myoblast',
           'LCL',
           'Bcell')

#rr <- rr[,n]

names(rr)[n]   <- namez
names(rr)[n+1] <- paste0('sim_',namez)

rForCC <- rr[,n,with=FALSE]
rOK    <- rr[,!grepl('^(expRT|simRT)',names(rr)), with=FALSE]

rr     <- rOK
r12    <- rr[rr$cs == cs2use,]

xMin <- 0
xMax <- ceiling(max(r12$to/1000000))

#### DONE WITH DATA IMPORT

#aGC  <-aggregate(r12$pcGC,   by=list(pos=round((r12$from+r12$to)/2/250000)*250000/1000000),FUN=mean)
aGC <- data.frame(pos=rollapply((r12$from+r12$to)/2,100,mean)/1000000,
                  x=rollapply(r12$pcGC,100,mean))

#aHiC <-aggregate(r12$hiCZyg, by=list(pos=round((r12$from+r12$to)/2/250000)*250000/1000000),FUN=mean)
aHiC <- data.frame(pos=rollapply((r12$from+r12$to)/2,100,mean)/1000000,
                  x=rollapply(r12$hiCZyg,100,mean))

aOri <-aggregate(r12$origins,by=list(pos=round((r12$from+r12$to)/2/250000)*250000/1000000),FUN=mean)

aGC$type <- 'GC'
aOri$type <- 'Ori'
aHiC$type <- 'HiC'

#aAll <- rbind(aGC,aHiC,aOri)

noXPos  <- theme(axis.text.x=element_blank())
noXLine <- theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())
noYPos  <- theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())
noMarg  <- theme(plot.margin = unit(c(-0.1,0,-0.1,0),'cm'))
bMarg   <- theme(plot.margin = unit(c(-0.1,0,0,0),'cm'))
smMarg  <- theme(plot.margin = unit(c(-0.1,0,0,0),'cm'))
noLeg   <- theme(legend.position='none')

#pMei <- plotRT(r12,'expRT_MeiS_ALL_S1_2to4C_SCP3yH2AX_NA_mm10',smoothSpan = 200) +
#pMei <- plotRT(r12,'expRT_MeiS_ALL_S5_2to4C_STRA8_DMRT1_mm10',smoothSpan = 100) +
pMei <- plotRT(r12,'MeiS(a)',smoothSpan = 100,yLabel = 'Mei-S') +
  noXPos +
  noMarg +
  xlab('')

pESC <- plotRT(r12,'ESC(a)',smoothSpan = 100,yLabel = 'ESC') +
  noXPos + noMarg  + xlab ('')

## Get Ori Plot
hmOri <- ggplot(aOri,aes(x=pos,y=1,fill=x^0.5)) +
  geom_tile(size=0) +
  coord_cartesian(expand = FALSE,xlim=c(xMin,xMax)) +
  noLeg +
  scale_fill_gradient(low=alpha('white',1),high='black') +
  noXPos +
  smMarg + noXLine +
  scale_y_continuous(breaks=c(1),labels=('Ori')) +
  xlab('') + ylab('') +
  theme(axis.ticks.y = element_blank())

## Prune GC outliers
aGC <- aGC[aGC$pos > 3,]
gcMed <- as.numeric(quantile(aGC$x,.5))
gcQ5  <- as.numeric(quantile(aGC$x,.05))
gcQ95 <- as.numeric(quantile(aGC$x,.95))

aGC$x[aGC$x>gcQ95] <- NA
aGC$x[aGC$x<gcQ5]  <- NA
aGC                <- aGC[!is.na(aGC$x),]
aGC$l2             <- log2(aGC$x/median(aGC$x))

aGC$up   <- aGC$l2
aGC$down <- aGC$l2

aGC$up[aGC$up <= 0]     <- 0
aGC$down[aGC$down >= 0] <- 0

#ylab(bquote(log[2]*'(GC/mean)')) +
#geom_smooth(lwd=.3,span=0.15,color='grey50') +#
#geom_line(aes(x=pos,y=l2),color='black',lwd=.3)+
lnGC <-ggplot(aGC,aes(x=pos,y=l2)) +
  geom_area(aes(x=pos,y=up),fill='red',lwd=.3,color=darkenColor('red',1.5)) +
  geom_area(aes(x=pos,y=down),fill='pink',lwd=.3,color=darkenColor('pink',1.5)) +

  coord_cartesian(expand = FALSE,xlim=c(xMin,xMax)) +
  geom_hline(yintercept=0,lwd=.2,lty='dashed')+
  noLeg + noXPos + noXLine + noMarg +
  ylab(bquote(atop("",atop(log[2],"(GC/med)")))) +
  xlab('') + 
  scale_y_continuous(breaks=c(-0.2,-0.1,0,0.1,0.2),labels=c('',-0.1,0,0.1,'')) +
  theme(axis.line.y=element_line(color='black',size=0.2),
        axis.title.y  = element_text(size = 10,
                                     margin = margin(t = 0, r = 0, b =0, l = -2, unit = "mm")))

aHiC$up   <- aHiC$x
aHiC$down <- aHiC$x

aHiC$up[aHiC$up <= 0]     <- 0
aHiC$down[aHiC$down >= 0] <- 0

hiCMin <- min(aHiC$x)*1.1
hiCMax <- max(aHiC$x)*1.1

#geom_line(aes(x=pos,y=x),color='black',lwd=.3)+
  
hmHiC <- ggplot(aHiC,aes(x=pos,y=1,fill=x)) +
  geom_area(aes(x=pos,y=up),fill='firebrick',alpha=.6,lwd=.3,color=darkenColor('firebrick',1.5)) +
  geom_area(aes(x=pos,y=down),fill='dodgerblue1',alpha=.6,lwd=.3,color=darkenColor('dodgerblue1',1.5))+
  geom_hline(yintercept=0,lwd=.1,lty='dashed') +
  coord_cartesian(expand = FALSE,xlim=c(xMin,xMax),ylim=c(hiCMin,hiCMax)) + noLeg +
  noXPos + bMarg +
  xlab(paste0('Position on chromosome ',csNum,' (Mb)')) +
  theme(axis.text.x=element_text(size=7)) +
  ylab("Hi-C\n(eig.)") +
  scale_y_continuous(breaks=c(-0.05,0,0.05),labels=c(-0.05,0,0.05)) +
  annotate(geom='text',x=45,y=-0.04,label='A',size=8*5/14,fontface='bold',color='firebrick') +
  annotate(geom='text',x=32,y=0.04,label='B',size=8*5/14,fontface='bold',color='dodgerblue2')

#theme(axis.ticks.x = element_blank(),
 #     axis.line.x = element_blank()) +

hmBlank <- ggplot(aHiC,aes(x=pos,y=1,fill=x)) + coord_cartesian(expand = FALSE,xlim=c(xMin,xMax)) + noLeg +
  scale_fill_gradient2(low='red',high='dodgerblue1',mid='black') + noYPos + smMarg + ylab('') +
  xlab(paste0('Position on chromosome ',csNum,' (Mb)'))

#hmHiC <- hmHiC + theme(axis.title.y=element_text(vjust=0.5,angle=0))
hmOri <- hmOri + theme(axis.title.y=element_text(vjust=0.5,angle=0))
#lnGC  <- lnGC + theme(axis.title.y=element_text(vjust=0.5 ,angle=0))

gAll <- ggarrange(pMei + noXLine,
                  hmOri,
                  pESC + noXLine,
                  lnGC + theme(plot.margin = margin(c(0,0,0,0),'cm')),
                  hmHiC,
          ncol=1,nrow=5,
          heights=c(4.2,1.1,4.2,2.8,2.8),
          labels = c('E','','F','H',''),
          font.label = list(size=8),
          hjust=0,vjust=1,
          align='v')

gSort <- plotSorting()

gX <- plotCCHeatMap(rForCC)

gABC <- ggarrange(gSort$gC,
                  gSort$gA + coord_cartesian(xlim=c(200,40000),
                                           ylim=c(20,70000),
                                           expand=FALSE),
                  gSort$gB,
                  ncol=3,nrow=1,
                  labels = c('A','B','C'),
                  font.label = list(size=8),
                  align='h',
                  hjust=0,vjust=1)

mypng <- readPNG(source = 'MeiS_purity_assessment_microscopy.png')
gD  <- rasterGrob(mypng, interpolate=TRUE)

gABCD <- ggarrange(gABC,
                   gD,
                  ncol=1,nrow=2,
                  heights=c(7,6),
                  labels = c('','D'),
                  font.label = list(size=8),
                  align='h',
                  hjust=0,vjust=1)

gEFGH <- ggarrange(gAll,gX,
                widths=c(9,6),
                ncol=2,nrow=1,
                labels = c('','G'),
                font.label = list(size=8),
                hjust=0,vjust=1)

gFig <- ggarrange(gABCD,gEFGH,ncol=1,nrow=2,heights=c(13,12))

scaleMe <- 1

ggsave('Pratto_et_al_Figure3.png',gFig,dpi=500,width=6.5*scaleMe,height=6.5*scaleMe)
ggsave('Pratto_et_al_Figure3.pdf',gFig,width=6.5*scaleMe,height=6.5*scaleMe)

cs2use <- 'chr12'
r12    <- rr[rr$cs == cs2use,]
pLst <- list()
### Associated supplement
for (nm in namez[2:length(namez)]){
  pLst[[nm]] <- plotRT(r12,nm,smoothSpan = 100,yLabel = nm) +
    noXPos +
    noMarg +
    annotate(geom='label',label=paste0(nm),x=-Inf,y=Inf,size=7*5/14,label.size=NA,hjust=0,vjust=1,) +
    ylab('') +
    xlab('')
}

xOK <- theme(axis.text.x=element_text(size=6),
             axis.ticks.x=element_line(),
             axis.title.x = element_text(size=8),
             plot.margin=unit(c(0,0.1,0,0),'cm'))

nLst <- length(pLst)
pLst[[nLst-1]] <- pLst[[nLst-1]] + xlab(paste0('\nPosition on ',cs2use,' (Mb)')) + xOK
pLst[[nLst]]   <- pLst[[nLst]]   + xlab(paste0('\nPosition on ',cs2use,' (Mb)')) + xOK

gALLRT <- grid.arrange(left=grid::textGrob(bquote('Replication timing; '*log[2]*'('*cov[n]*'/median)'),gp=gpar(fontsize=8),rot=90,vjust=1),
                       ggarrange(plotlist = pLst,align='hv',ncol=2,nrow=ceiling(length(pLst)/2)))

scaleMe <- 0.9
ggsave('Pratto_et_al_SuppFig_RT_all_celltypes.png',gALLRT,dpi=400,width=7*scaleMe,height=8*scaleMe)
ggsave('Pratto_et_al_SuppFig_RT_all_celltypes.pdf',gALLRT,width=7*scaleMe,height=9*scaleMe)
