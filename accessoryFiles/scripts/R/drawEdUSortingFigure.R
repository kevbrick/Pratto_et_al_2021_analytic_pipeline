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

getSort1Data <- function(fcsIn = '2019-08-16_Edu-AF488_001_006.fcs', limz = NULL) {
  dSort        <- importFCS(fcsIn)
  names(dSort)[names(dSort) == 'BV421.A'] <- 'DAPI'
  names(dSort)[names(dSort) == 'G575.A'] <- 'STRA8'
  names(dSort)[names(dSort) == 'B530.A'] <- 'EdU'
  names(dSort)[names(dSort) == 'R660.A'] <- 'DMRT1'
  
  dSort$DAPI1C  <- dSort$DAPI  > 35000 & dSort$DAPI < 60000
  dSort$DAPI2C  <- dSort$DAPI  > 70000 & dSort$DAPI < 90000
  dSort$DAPIRep  <- dSort$DAPI  > 110000 & dSort$DAPI < 150000
  dSort$DAPI4C  <- dSort$DAPI  > 170000 & dSort$DAPI < 220000
  dSort$DAPIok  <- dSort$DAPI  > 30000 & dSort$DAPI < 220000 & dSort$BV421.W < 100000
  dSort$fok <- dSort$FSC.A > 40000 & dSort$FSC.H < 100000
  dSort$STRA8ok  <- dSort$STRA8  > 7000
  dSort$DMRT1ok <- dSort$DMRT1 > 800
  dSort$EdUp <- dSort$EdU > 500
  
  
  dSort$DSok    <- dSort$DAPIok 
  
  dSort$allOK             <- 'Other'
  dSort$allOK[dSort$DSok] <- 'DAPI_Gated'
  dSort$ploidy <- 'undet'
  dSort$ploidy[dSort$DAPI1C] <- '1C'
  dSort$ploidy[dSort$DAPI2C] <- '2C'
  dSort$ploidy[dSort$DAPIRep] <- 'Rep'
  dSort$ploidy[dSort$DAPI4C] <- '4C'
  
  dSortOK <- dSort[dSort$DSok & dSort$fok,]
  
  return(list(all = dSort, ok=dSortOK))
}

dfSortData <- getSort1Data()

dfRepStatus <- dfSortData$ok[dfSortData$ok$DAPIRep,]

dfRepStatus$stage <- "Unknown"
dfRepStatus$stage[dfRepStatus$DMRT1 > 3000 & dfRepStatus$DMRT1 < 15000 & dfRepStatus$STRA8 < 7000  & dfRepStatus$STRA8 > 1000 ]  <- 'SpgU'
dfRepStatus$stage[dfRepStatus$DMRT1 > 15000 & dfRepStatus$DMRT1 < 40000 & dfRepStatus$STRA8 < 7000  & dfRepStatus$STRA8 > 1000 ] <- 'SpgI'
dfRepStatus$stage[dfRepStatus$DMRT1 > 3000 & dfRepStatus$DMRT1 < 40000 & dfRepStatus$STRA8 < 50000  & dfRepStatus$STRA8 > 10000] <- 'SpgB'
dfRepStatus$stage[dfRepStatus$DMRT1 < 2500 & dfRepStatus$DMRT1 > 500 & dfRepStatus$STRA8 < 50000  & dfRepStatus$STRA8 > 10000]   <- 'MeiS'

colz   <- c("SpgB" = "#009E73", 
            "SpgU"="#CC79A7", 
            "SpgI"="#D55E00", 
            "MeiS" = "#444444",
            "Unknown" = "#cccccc", 
            "Others" = "#cccccc", 
            "Other2" = "#56B4E9", 
            "Other3" = "#0072B2",
            "Other4" = "#F0E442")

## Panel A : Ploidy ---------------------------------------------------------------------
dfPloidySubset <- dfSortData$ok

gPloidy <- ggplot(dfPloidySubset,
       aes(x=DAPI/1000,
           y=EdU)) + 
  geom_point(alpha=.05,size=.1,color='grey10') + 
  scale_y_log10(breaks=c(10,100,1000,10000),labels=fancy_scientific) + 
  annotation_logticks(sides='l',
                      size=0.2,
                      short=unit(0.050,'cm'),
                      mid=unit(0.075,'cm'),
                      long=unit(0.100,'cm'))+
  geom_hline(yintercept=500, color='red', lwd=.3, lty='dashed')+
  geom_text(data=data.frame(DAPI=120000,
                            EdU=3000),
            label='EdU positive\n1.8%',
            vjust=1,
            hjust=0.5,
            check_overlap = TRUE,
            size=7*5/14,
            color='red')+
  geom_text(data=data.frame(DAPI=50000,
                            EdU=30),
            label='1C',
            vjust=1,
            hjust=0.5,
            check_overlap = TRUE,
            size=7*5/14,
            color='black')+
  geom_text(data=data.frame(DAPI=100000,
                            EdU=100),
            label='2C',
            vjust=1,
            hjust=0.5,
            check_overlap = TRUE,
            size=7*5/14,
            color='black')+
  geom_text(data=data.frame(DAPI=170000,
                            EdU=100),
            label='4C',
            vjust=1,
            hjust=0.5,
            check_overlap = TRUE,
            size=7*5/14,
            color='black')+
  xlab(expression('DAPI fluorescence (x10'^3*')')) +
  ylab(expression('EdU fluorescence')) + 
  theme(legend.position='none')+ coord_cartesian(xlim = c(20,220),ylim=c(20,4000))

dfDensity <- dfPloidySubset %>% add_tally(name='tot') %>% 
  mutate(DAPIbin=round(DAPI/1000)) %>%
  group_by(tot,DAPIbin) %>% 
  count() %>% mutate(pc=n/tot*100)

gDAPIdensity <- ggplot(dfDensity,aes(x=DAPIbin,y=pc)) + 
  scale_y_log10() + 
  coord_cartesian(ylim=c(0.01,4)) +
  ylab(expression('Density')) +
  annotation_logticks(sides='l',
                      size=0.2,
                      short=unit(0.050,'cm'),
                      mid=unit(0.075,'cm'),
                      long=unit(0.100,'cm')) +
  geom_smooth(span=.3,color='black',lwd=.3) + 
  xlab(expression('DAPI fluorescence (x10'^3*')')) +
  ylab('Nuclei (%)')

gPloidyWithDensity <- ggarrange(gPloidy,gDAPIdensity,
                                ncol=1,
                                heights=c(4,1),
                                labels=c('A',''),
                                font.label = list(size=8,fontface='bold'),
                                hjust=0,vjust=1)

## Panel B: Sort -------------------------------------------------------------------------

dfLabels <- dfRepStatus %>% group_by(stage) %>% 
  add_tally(name = 'tot') %>% 
  group_by(stage,EdUp,tot) %>% 
  tally() %>% 
  dplyr:::filter(EdUp) %>% 
  mutate(pc=n/tot*100,
         lbl=paste0(stage,"\n",round(pc),"% EdU positive"))

gSortingSchema <- ggplot(dfRepStatus,
                         aes(x=STRA8, y=DMRT1, color=EdUp)) + 
  geom_point(size=.3,alpha=0.2,lwd=.2) +
  scale_x_log10(breaks=c(1000,10000,100000), labels=fancy_scientific) + 
  scale_y_log10(breaks=c(100,1000,10000,100000), labels=fancy_scientific) + 
  annotate(geom='rect',
           xmax=7000,xmin=1000,
           ymax=40000,ymin=15000,
           lwd=.4, lty='dashed',
           fill=NA,
           color=colz[names(colz)=='SpgI'])+
  annotate(geom='rect',
           xmax=7000,xmin=1000,
           ymax=15000,ymin=3000,
           lwd=.4, lty='dashed',
           fill=NA,
           color=colz[names(colz)=='SpgU'])+
  annotate(geom='rect',
           xmax=50000,xmin=10000,
           ymax=40000,ymin=3000,
           lwd=.4, lty='dashed',
           fill=NA,
           color=colz[names(colz)=='SpgB'])+
  annotate(geom='rect',
           xmax=50000,xmin=10000,
           ymax=2500,ymin=500,
           lwd=.4, lty='dashed',
           fill=NA,
           color=colz[names(colz)=='MeiS'])+
  annotation_logticks(sides='bl',
                      size=0.2,
                      short=unit(0.050,'cm'),
                      mid=unit(0.075,'cm'),
                      long=unit(0.100,'cm'))+ 
  scale_color_manual(values=c('grey70','red')) + 
  theme(legend.title=element_blank(),
        legend.position='none',
        legend.justification=c(1,0),
        legend.direction = 'horizontal',
        legend.key.size = unit(.2,'cm'),
        legend.background=element_blank()) + 
  xlab('STRA8 fluorescence') +
  ylab('DMRT1 fluorescence') +
  coord_cartesian(ylim=c(80,100000),xlim=c(500,60000))+
  annotate(geom='text',
           x=50000,y=500,
           label=dfLabels$lbl[dfLabels$stage == 'MeiS'],
           hjust=1,vjust=1.5,
           size=7*5/14,
           color=colz[names(colz)=='MeiS']) +
  annotate(geom='text',
           x=1000,y=3000,
           label=dfLabels$lbl[dfLabels$stage == 'SpgU'],
           hjust=0,vjust=1.5,
           size=7*5/14,
           color=colz[names(colz)=='SpgU'])+ 
  annotate(geom='text',
           x=50000,y=40000,
           label=dfLabels$lbl[dfLabels$stage == 'SpgB'],
           hjust=1,vjust=-.5,
           size=7*5/14,
           color=colz[names(colz)=='SpgB']) +
  annotate(geom='text',
           x=1000,y=40000,
           label=dfLabels$lbl[dfLabels$stage == 'SpgI'],
           hjust=0,vjust=-0.5,
           size=7*5/14,
           color=colz[names(colz)=='SpgI'])

## Panel C: Replicating cell proportions ---------------------------------------------
dfReplicatingCells <- dfRepStatus %>% dplyr::count(stage,EdUp) %>% 
  dplyr::filter(EdUp=="TRUE")   %>%
  dplyr::mutate(allrep=sum(n))  %>%  
  dplyr::mutate(percentage=round(n/allrep*100,1),
                label=paste0(stage,"\n",percentage," %"),
                pcLabel=paste0(percentage," %"),
                Stage=ifelse(stage=='Unknown','Others',stage)) %>%
  dplyr::select(Stage, percentage, label, pcLabel)

dfReplicatingCells$Stage <- factor(dfReplicatingCells$Stage, 
                                   levels = c("MeiS","SpgU","SpgB","SpgI","Others"))

dfReplicatingCells$ypos <- c(11,75,17.5,26,2)

gBar <- ggplot(dfReplicatingCells, 
               aes(x=Stage,
                   y=percentage,
                   fill=Stage)) +
  geom_bar(stat="identity", width=1, color="black",lwd=.2) +
  theme(legend.position="none") +
  scale_fill_manual(values=colz) + 
  geom_text(aes(y = percentage+5, 
                x=Stage, 
                label = pcLabel), 
            color = "black", 
            size=7*5/14,
            vjust=1) +
  xlab('') + 
  coord_cartesian(ylim=c(0,75),xlim=c(0.45,5.55),expand=FALSE) +
  ylab('Replicating cells in testis (%)')

gPie <- ggplot(dfReplicatingCells, aes(x="",
                                       y=percentage,
                                       fill=Stage)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_label_repel(aes(x=1,label = label), color='black',position=position_stack(), size=7*5/14) +
  scale_fill_manual(values=colz) + 
  scale_y_reverse()

## Merge ----------------------------------------------------------------------------------
gBC <- ggarrange(gSortingSchema,gBar,
                 ncol=1,nrow=2,
                 labels=c('B','C'),
                 heights=c(2,1.3),
                 font.label = list(size=8,fontface='bold'),
                 hjust=0,vjust=1)

gABC <- ggarrange(gPloidyWithDensity,gBC,
                  ncol=2,nrow=1,
                  widths=c(5,4))

scale <- .9

ggsave('Pratto_et_al_SuppFig_TestisCellPops.png',gABC,width=6.5*scale,height=5*scale)
ggsave('Pratto_et_al_SuppFig_TestisCellPops.pdf',gABC,width=6.5*scale,height=5*scale)

# ## Merge ----------------------------------------------------------------------------------
# gBC <- ggarrange(gSortingSchema,gPie,
#                  ncol=1,nrow=2,
#                  labels=c('B','C'),
#                  heights=c(2,1.3),
#                  font.label = list(size=8,fontface='bold'),
#                  hjust=0,vjust=1)
# 
# gABC <- ggarrange(gPloidyWithDensity,gBC,
#                   ncol=2,nrow=1,
#                   widths=c(5,4))
# 
# scale <- 1
# 
# ggsave('Pratto_et_al_SuppFig_TestisCellPops_Pie.png',gABC,width=6.5*scale,height=5*scale)
# ggsave('Pratto_et_al_SuppFig_TestisCellPops_Pie.pdf',gABC,width=6.5*scale,height=5*scale)