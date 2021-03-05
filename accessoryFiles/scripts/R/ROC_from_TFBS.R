## Load custom functions & H3 data
source('accessoryFiles/scripts/R/genericFunctions.R')
source('accessoryFiles/scripts/R/repFunctions.R')

theme7point()

nDecimals <- 1

calc_auc <- function(X, Y){
  # inputs already sorted, best scores first
  dY <- c(diff(Y), 0)
  dX <- c(diff(X), 0)
  sum(Y * dX) + sum(dY * dX)/2
}

m <- fread('motifs.bed')

names(m) <- c('cs','from','to','name','score','strand','peak')
m$atPeak <- m$peak>0
m$sc <- round(m$score,nDecimals)

mDF <- as.data.frame(m %>% group_by(sc) %>% summarise(nMotifs = sum(atPeak>0),FPmotifs = sum(atPeak == 0)))

#ggplot(mDF,aes(x=sc,y=TP/(TP+FP))) + geom_point()

p <- fread('peaks_with_motifs.bed')
names(p) <- c('cs','from','to','score')

p$score <- as.numeric(p$score)
p <- p[!is.na(p$score),]

p$sc <- round(p$score,nDecimals)

pDF <- as.data.frame(p %>% group_by(sc) %>% summarise(TP = sum(as.numeric(score)>-999)))

pDF$detected <- rev(cumsum(rev(pDF$TP)))
pDFU <- pDF[pDF$detected <= 12000,]
pDFU <- pDFU[!pDFU$sc == -999, ]

pX2 <- plyr:::join(pDFU,mDF,by="sc",type='left')

pX2$FN <- max(pX2$detected) - pX2$detected

pX2$sens <- pX2$detected / (pX2$detected+pX2$FN)

pX2$cumTPmotifs <- rev(cumsum(rev(pX2$nMotifs)))
pX2$cumFPmotifs <- rev(cumsum(rev(pX2$FPmotifs)))
pX2$TN <- max(pX2$cumFPmotifs) - pX2$cumFPmotifs
pX2$spec <- pX2$TN/(pX2$TN+pX2$cumFPmotifs)

pX2$name <- gsub(pattern = ';.+$',replacement = '',m$name[1])
pX2$AUC  <- round(calc_auc(1-pX2$sens,pX2$spec),2)

plotDF <- pX2[,c('name','sc','spec','sens','AUC')]

gROC <- ggplot(plotDF,aes(x=1-spec,y=sens)) +
  geom_line(lwd=.3,color='red') +
  geom_point(size=.4,color='darkred') +
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color='grey60',lwd=.2) +
  coord_cartesian(xlim=c(-0.02,1.025),ylim=c(-0.02,1.025),expand=FALSE) +
  xlab('1 - Specificity') +
  ylab('Sensitivity') +
  annotate(geom='text',label=paste0(plotDF$name[1],"\nAUC = ",plotDF$AUC[1]),x=.75,y=.25,size=8*5/14)

ggsave(paste0(plotDF$name[1],'_ROC.png'),gROC,height=2,width=2)
ggsave(paste0(plotDF$name[1],'_ROC.pdf'),gROC,height=2,width=2)

write.table(x = plotDF, file = paste0(plotDF$name[1],'_ROC.tab'), quote = FALSE, row.names = FALSE, col.names = TRUE)
