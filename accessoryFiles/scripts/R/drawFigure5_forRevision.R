imgdir  <- '.'

## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')
source('accessoryFiles/scripts/R/simulateReplicationTiming_allCS.R')
source('accessoryFiles/scripts/R/repliSim_loadModules.R')

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(data.table)

theme7point()

tInit      <- fread('rep_v_rec_MM10.tab',header=TRUE)
tColz      <- names(tInit)[grep('exp',names(tInit))]
tS         <- apply(tInit[,..tColz],1,function(x){sum(abs(x) < 0.001)>1})
t          <- tInit[!tS,]

randomReps <- 0
pointSzie2Use <- 0.2
nQuantiles  <- 250
nSmoothSpan <- 0.8
###############################################
nuDHH      <- plotFeaturesVRTMulti(iDF     = t,
                                   iFields = c("h3k4m3hh","SSDShh"),
                                   names   = c('Pre-DSB mark \n(h3k4m3) ',
                                               'DSB + repair \n(Dmc1 SSDS) '),
                                   iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                   q       = nQuantiles,
                                   rev     = TRUE,
                                   smoothSpan=nSmoothSpan,
                                   pointSz = pointSzie2Use,
                                   nRandom = randomReps,
                                   pCol    = 'darkorange',
                                   ylbl    = bquote(B6^"h/h"*" "))

a  <- annotate(geom='text',x=Inf,y=Inf,label=' Pre-DSB mark ',color='orange',hjust=1,vjust=1,size=7*5/14)
a2 <- annotate(geom='text',x=Inf,y=Inf,label='\n DSB & repair ',color='forestgreen',hjust=1,vjust=1,size=7*5/14)
gHH <- nuDHH$figAlt + coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) + a + a2

nuDB6         <- plotFeaturesVRTMulti(iDF     = t,
                                    iFields = c("h3k4m312DPPR1","spo11","SSDST1"),
                                    names   = c('Pre-DSB mark \n(h3k4m3) ',
                                                'DSB repair (spo11)',
                                                'DSB + repair \n(Dmc1 SSDS) '),
                                    iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                    q       = nQuantiles,
                                    rev     = TRUE,
                                    smoothSpan=nSmoothSpan,
                                    pointSz = pointSzie2Use,
                                    nRandom = randomReps,
                                    pCol    = 'darkorange',
                                    ylbl    = bquote(B6^"wt/wt"*" "))

a  <- annotate(geom='text',x=Inf,y=Inf,label=' Pre-DSB mark ',color='orange',hjust=1,vjust=1,size=7*5/14)
a2 <- annotate(geom='text',x=Inf,y=Inf,label='\n DSBs ',color='magenta',hjust=1,vjust=1,size=7*5/14)
a3 <- annotate(geom='text',x=Inf,y=Inf,label='\n\n DSB & repair ',color='forestgreen',hjust=1,vjust=1,size=7*5/14)
gB6 <- nuDB6$figAlt + coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) + a + a2 + a3

#############################################
t$myersCONCOhh <- t$myersCOHum + t$myersNCOHUM
nuCONCO  <- plotFeaturesVRTMulti(iDF     = t,
                                 iFields = c("h3k4m3hh","myersCONCOhh"),                                  
                                 names   = c('Pre-DSB mark (h3k4m3) ','IH repair'),
                                   iRT     = "simRT_MeiS_BxCMID_R1_2to4C_SCP3_NA_mm10",
                                   q       = nQuantiles,
                                   rev     = TRUE,
                                   smoothSpan=nSmoothSpan,
                                   pointSz = pointSzie2Use,
                                   nRandom = randomReps,
                                   pCol    = 'darkorange',
                                   ylbl    = bquote(CASTxB6^"h/h"*"\n"))

a  <- annotate(geom='text',x=Inf,y=Inf,label=' Pre-DSB mark ',color='orange',hjust=1,vjust=1,size=7*5/14)
a2 <- annotate(geom='text',x=Inf,y=Inf,label='\n Inter-homolog repair' ,color='black'  ,hjust=1,vjust=1,size=7*5/14)

gCO <- nuCONCO$figAlt + coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) + a + a2 

#############################################
nuAffy   <- plotFeaturesVRTMulti(iDF     = t,
                                 iFields = c("affySeq"),#,"Prdm9ChIPSeq"),
                                 names   = c('affySeq'),#'PRDM9ChIPSeq'),
                                 iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                 q       = nQuantiles,
                                 smoothSpan=nSmoothSpan,
                                 pointSz = pointSzie2Use,
                                 rev     = TRUE,
                                 nRandom = randomReps,
                                 pCol    = 'darkorange',
                                 ylbl    = bquote(B6^"wt/wt"*" "))

#a  <- annotate(geom='text',x=Inf,y=Inf,label=' PRDM9 in vitro',color='royalblue1',hjust=1,vjust=1,size=7*5/14)
a2  <- annotate(geom='text',x=Inf,y=Inf,label=' PRDM9 in vitro',color='dodgerblue4',hjust=1,vjust=1,size=7*5/14)

gAffy <- nuAffy$figAlt + coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) + a2

#
#############################################
nuPRKO   <- plotFeaturesVRTMulti(iDF     = t,
                                 iFields = c("SSDSPrKO"),
                                 names   = c("SSDSPrKO"),
                                 iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                 q       = nQuantiles,
                                 smoothSpan=nSmoothSpan,
                                 pointSz = pointSzie2Use,
                                 rev     = TRUE,
                                 nRandom = randomReps,
                                 pCol    = 'darkorange',
                                 ylbl    = bquote("B6 "*Prdm9^"-/-"*" "))

#a  <- annotate(geom='text',x=Inf,y=Inf,label='DSB rate (Spo11)',color='royalblue1',hjust=-0.1,vjust=1,size=7*5/14)

nm <- 'DSB + repair\n (Prdm9 KO)'
a2  <- annotate(geom='text',x=Inf,y=Inf,label=nm,color='forestgreen',hjust=1,vjust=1,size=7*5/14)

gPRKO <- nuPRKO$figAlt + coord_cartesian(ylim=c(-2,2)) + 
  scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),
                     labels=c('Early','Mid','Late')) + 
  a2

#############################################
nuHOP2KO   <- plotFeaturesVRTMulti(iDF     = t,
                                 iFields = c("SSDShop2"),
                                 names   = c("SSDShop2"),
                                 iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                 q       = nQuantiles,
                                 smoothSpan=nSmoothSpan,
                                 pointSz = pointSzie2Use,
                                 rev     = TRUE,
                                 nRandom = randomReps,
                                 pCol    = 'darkorange',
                                 ylbl    = bquote("B6 "*Hop2^"-/-"*" "))

#a  <- annotate(geom='text',x=Inf,y=Inf,label='DSB rate (Spo11)',color='royalblue1',hjust=-0.1,vjust=1,size=7*5/14)

nm <- 'DSB + repair\n (Hop2 KO)'
a2  <- annotate(geom='text',x=Inf,y=Inf,label=nm,color='forestgreen',hjust=1,vjust=1,size=7*5/14)

gHOP2KO <- nuHOP2KO$figAlt + coord_cartesian(ylim=c(-2,2)) + 
  scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),
                     labels=c('Early','Mid','Late')) + 
  a2

#############################################
nuRPA     <- plotFeaturesVRTMulti(iDF     = t,
                                  iFields = c("RPA","SSDST1","SSDST2","SSDShop2","SSDSPrKO","spo11"),
                                  names   = c("RPA","SSDST1","SSDST2","SSDShop2","SSDSPrKO","spo11"),
                                  iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
                                  q       = nQuantiles,
                                  smoothSpan=nSmoothSpan,
                                  pointSz = pointSzie2Use,
                                  rev     = TRUE,
                                  nRandom = randomReps,
                                  pCol    = 'darkorange',
                                  ylbl    = "B6 ")

nuCONCOSupp     <- plotFeaturesVRTMulti(iDF     = t,
                                  iFields = c('myersCOHum','myersNCOHUM'),
                                  names   = c('CO','NCO'),
                                  iRT     = "simRT_MeiS_BxCMID_R1_2to4C_SCP3_NA_mm10",
                                  q       = nQuantiles,
                                  smoothSpan=nSmoothSpan,
                                  pointSz = pointSzie2Use,
                                  rev     = TRUE,
                                  nRandom = randomReps,
                                  pCol    = 'darkorange',
                                  ylbl    = bquote(B6xCAST^"h/h "))

nSlope <- 1.8

gRPA <- nuRPA$fig +
  coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) +
  theme(panel.grid=element_line(size=.25,linetype='dotted',color='grey70')) +
  geom_line(data=data.frame(x=c(0,nQuantiles),y=c(nSlope,-nSlope)),aes(x=x,y=y,fill=NA),color='grey60') +
  geom_hline(yintercept=0,lwd=.3,alpha=.3) + facet_wrap(~name,nrow=2)

gCNC <- nuCONCOSupp$fig +
  coord_cartesian(ylim=c(-2,2)) + 
  scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late')) +
  theme(panel.grid=element_line(size=.25,linetype='dotted',color='grey70')) +
  geom_line(data=data.frame(x=c(0,nQuantiles),y=c(nSlope,-nSlope)),aes(x=x,y=y,fill=NA),color='grey60') +
  geom_hline(yintercept=0,lwd=.3,alpha=.3) + facet_wrap(~name,nrow=2)

gS <- ggarrange(gRPA,gCNC,widths=c(3,1),ncol=2,
                align='hv',
                labels=c('F',''),
                hjust = 0, vjust = 1,
                font.label = list(size=8))
# N=1.4
# ggsave(plot = gS,filename = 'Pratto_figure4_Prdm9KOandRPASSDSSupplement.png',width=4.7*N,height=2.5*N)
# ggsave(plot = gS,filename = 'Pratto_figure4_Prdm9KOandRPASSDSSupplement.pdf',width=4.7*N,height=2.5*N)

################################################
gBCD <- ggarrange(gB6,gAffy,gHH,gCO,
          align='hv',
          labels=c('B','C','D','E'),
          hjust = 0, vjust = 1,
          font.label = list(size=8),
          nrow=2,ncol=2,
          heights=c(5,5))

gBCDEFG <- ggarrange(gB6,gAffy,gHH,gHOP2KO,gPRKO,gCO,
                  align='hv',
                  labels=c('B','C','D','E','F','G'),
                  hjust = 0, vjust = 1,
                  font.label = list(size=8),
                  nrow=2,ncol=3,
                  heights=c(5,5))

gBCDE <- ggarrange(gB6,gAffy,gHH,gCO,
                     align='hv',
                     labels=c('B','C','D','E'),
                     hjust = 0, vjust = 1,
                     font.label = list(size=8),
                     nrow=2,ncol=2,
                     heights=c(5,5))

mypng <- readPNG(source = 'DSBschemaBIG.png')
gPic  <- rasterGrob(mypng, interpolate=TRUE)


gTop <- ggarrange(grid.arrange(gPic),
                  gBCDE,
                  widths=c(1.4,4),
                  labels = c('A',''),
                  hjust = 0, vjust = 1,
                  font.label = list(size=8),
                  ncol=2,nrow=1)

N=1

gFigure <- ggarrange(gTop,ggplot()+theme_void(),gS,ncol=1,nrow=3,heights=c(4,.2,3))

ggsave(plot = gFigure,filename = 'Pratto_et_al_figure5.png',width=7*N,height=7*N)
ggsave(plot = gFigure,filename = 'Pratto_et_al_figure5.pdf',width=7*N,height=7*N)

# gE <- ggarrange(gPRKO,
#                 gRPA,
#                 widths=c(1,1),
#                 labels = c('A','B'),
#                 hjust = 0, vjust = 1,
#                 font.label = list(size=8),
#                 ncol=2,nrow=1)

# gSlice <- drawRTvRepSlice(slice='chr12')
#
# gALL <- ggarrange(gSlice,gTop,nrow=2,ncol=1,
#                   labels = c('a',''),
#                   hjust = 0, vjust = 1,
#                   font.label = list(size=8),
#                   heights=c(1,4))
#
#
# N=0.9; ggsave(plot = gALL,filename = 'Pratto_figure4.png',width=7*N,height=3.5/3*4*N)
# N=0.9; ggsave(plot = gALL,filename = 'Pratto_figure4.pdf',width=7*N,height=3.5/3*4*N)

# #################################################
# 
# #############################################
# nuOth     <- plotFeaturesVRTMulti(iDF     = t,
#                                  iFields = c("pcGC","hiCZyg","h3k9m2"),#,"Prdm9ChIPSeq"),
#                                  names   = c("GC(%)","HiC","H3K9m2"),#'PRDM9ChIPSeq'),
#                                  iRT     = "simRT_MeiS_VERYEARLY_S1_2to4C_SCP3_NA_mm10",
#                                  q       = nQuantiles,
#                                  smoothSpan=.8,
#                                  pointSz = pointSzie2Use,
#                                  rev     = TRUE,
#                                  nRandom = randomReps,
#                                  pCol    = 'darkorange',
#                                  ylbl    = bquote(B6^"wt/wt"*" "))
# 
# gOth <- nuOth$fig + coord_cartesian(ylim=c(-2,2)) + scale_x_continuous(breaks = c(nQuantiles*.05,nQuantiles*.5,nQuantiles*.95),labels=c('Early','Mid','Late'))
# 
# #N=0.9; ggsave(plot = gOth,filename = 'Pratto_et_al_figure4_GCHiCHet.png',width=6*N,height=1.8/3*4*N)

#a  <- annotate(geom='text',x=Inf,y=Inf,label=' PRDM9 in vitro',color='royalblue1',hjust=1,vjust=1,size=7*5/14)
#a2  <- annotate(geom='text',x=Inf,y=Inf,label=' PRDM9 in vitro',color='dodgerblue4',hjust=1,vjust=1,size=7*5/14)
