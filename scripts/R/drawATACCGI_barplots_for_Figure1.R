## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

imgdir <- '.'

library(png)
library(ggpubr)
library(plyr)
library(ggplot2)
library(tidyverse)

theme7point()

## Panel E:
dfE <- read.table('figure1eData.tab',header=TRUE) 

dfPlotE <- dfE %>% mutate(pc=N/tot*100,
                         lbl=paste0(" ",round(N/tot*100,0),'%'))

names(dfPlotE)[3] <- 'ATAC\nSeq'

dfPlotE$name <- factor(dfPlotE$name, 
                      levels=rev(c('ATAC-CGI','ATACnoCGI','CGInoATAC','Neither')))

gBarsE <- ggplot(dfPlotE,aes(x=pc,y=name)) + 
  geom_bar(stat='identity',width = .6, fill='grey50',color='black',lwd=.3) + 
  xlab('Origins (%)') + 
  coord_cartesian(xlim=c(0,75),ylim=c(0.5,4.5),expand=FALSE) +
  theme(legend.position='none',
        axis.text.y=element_blank(),
        plot.margin=unit(c(0,0,0,0),'cm')) + 
  geom_text(aes(label=lbl),size=7*5/14,hjust=-0.5,color='black')+
  ylab('')

dfKeyE <- reshape2:::melt.data.frame(dfPlotE,id.vars = c('name','pc'),
                                    measure.vars = c('ATAC\nSeq','CGI'))

gKeyE <- ggplot(dfKeyE,aes(x=variable,y=name,fill=value)) + 
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

gE <- ggarrange(gKeyE,
                gBarsE,
                ncol=2,nrow=1,
                widths=c(1,4),
                align='h')

## Figure 1F:
dfF <- read.table('figure1fData.tab',header=TRUE,sep = "\t") 

dfPlotF <- dfF %>% mutate(pc=N/tot*100,
                         lbl=paste0("  ",round(N/tot*100,0),'%'))

dfPlotF$type <- factor(dfPlotF$type, 
                       levels=rev(c('ATAC & CGI','CGI no ATAC','ATAC no CGI')))

gBarsF <- ggplot(dfPlotF,aes(x=pc,y=type)) + 
  geom_bar(stat='identity',width = .5, fill='grey50',color='black',lwd=.3) + 
  xlab('At origins (%)') + 
  coord_cartesian(xlim=c(0,100),ylim=c(0.5,3.5),expand=FALSE) +
  theme(legend.position='none',
        plot.margin=unit(c(0,0,0,0),'cm')) + 
  geom_text(aes(label=lbl),size=7*5/14,hjust=-0.5,color='black')+
  ylab('')

gF <- ggarrange(ggplot() + theme_void(),
                gBarsF,
                ncol=2,nrow=1,
                widths=c(1,7),
                align='h')

ggarrange(gE,gBarsF,
          ncol=1,nrow=2,
          heights=c(5,5),
          labels = c('E','F'),
          vjust = 1,hjust=0,
          font.label = list(size=8,face='bold',vjust=1))

ggsave('Pratto_et_al_Figure1EF.png',width=1.9,height=2.4)
ggsave('Pratto_et_al_Figure1EF.pdf',width=1.9,height=2.4)
