library(mixtools)
library(seqinr)
library(ape)
library(edgeR)
library(DESeq)
library(lattice)
library(MASS)
library(gplots)
library(maptools)
library(dplyr)
library(tidyverse)

#########OR TREE PLOTS###################
plot.OR.tree <- function(OR.tree,edgecol=MyGray2,edgewid=0.5) { 
  plot(OR.tree,cex=1,type="unrooted",edge.color=edgecol,pch=21,show.tip.label=F,tip.color="black",edge.width=edgewid,no.margin=T)#,y.lim=c(-0.41,0.29)) x.lim=c(-0.2,0.4),y.lim=c(0.4,-0.3
  tiplabels(pch=21,cex=0.9,col=MyTransparent,bg=MyTransparent,lwd=0.7)
}

add.OR.tree.label <- function(OR.tree,label,pch=21,pt.col="azure2",bg.col="azure1",lwd=0.7,cex=1,frame="none") {
  length <- length(OR.tree$tip.label)
  pt.color <- rep(MyTransparent,length)
  bg.color <- rep(MyTransparent,length)
  tip.color <- rep(MyTransparent,length)
  for (i in 1:length) {
    if (OR.tree$tip.label[i] %in% label)
    { pt.color[i]<-pt.col;bg.color[i]<-bg.col;tip.color[i]<-bg.col }
  }   
  tiplabels(pch=pch,cex=cex,col=pt.color,bg=bg.color,lwd=lwd)
  #  tiplabels(OR.tree$tip.label,cex=cex,col=tip.color,bg=MyTransparent,frame="none")
}

auto.OR.tree <- function(result.odor.OR,odor) {
  to.edgeA <- subset(result.odor.OR,FDR<0.05&logFC>0)$Gene_name
  to.edgeB <- subset(result.odor.OR,FDR<0.001&logFC>0)$Gene_name
  
  edgeOdorA <- which.edge(OR.tree,to.edgeA)
  edgeOdorB <- which.edge(OR.tree,to.edgeB)
  edge.clr <- rep("gray87",length(OR.tree$edge)) #baisc lines "gray87"
  edge.clr[edgeOdorA] <- "blue"
  edge.clr[edgeOdorB] <- "red"
  edge.wdth <- rep(0.5,length(OR.tree$edge)) #gray line thickness 0.5
  edge.wdth[edgeOdorA] <- 3 #color line thickness 3
  edge.wdth[edgeOdorB] <- 3 #color line thickness 3
  
  plot.OR.tree(OR.tree,edge.clr,edge.wdth)
  add.OR.tree.label(OR.tree,result.odor.OR$Gene_name,21,"seashell3","ivory2",0.5) 
  #All ORs will be grey circles
  FDR.OR.A <- subset(result.odor.OR,FDR<0.05&logFC>0)$Gene_name
  FDR.OR.B <- subset(result.odor.OR,FDR<0.001&logFC>0)$Gene_name
  
  add.OR.tree.label(OR.tree,FDR.OR.A,21,"gray31","blue",0.5,1)
  add.OR.tree.label(OR.tree,FDR.OR.B,21,"gray31","red",0.5,1)
  
  plot.new();legend("topleft",c(paste(odor,"ORs p<0.001"),paste(odor,"ORs p<0.05")), pch=21,col="gray31",pt.bg=c("red","blue"),bty="n",pt.lwd=0.5,pt.cex=1,cex=1)  
}


pdf("STAR_ACT_TMT.OR.Tree.pdf",width=7, height=7)
auto.OR.tree(result.Ser_0.01pACT.STAR.OR, "0.01% Acetophenone")
auto.OR.tree(result.Ser_0.1pACT.STAR.OR, "0.1% Acetophenone")
auto.OR.tree(result.Ser_1pACT.STAR.OR,"1% Acetophenone")
auto.OR.tree(result.Ser_100pACT.STAR.OR,"100% Acetophenone") 
auto.OR.tree(result.Ser_0.01pTMT.STAR.OR,"0.01% TMT")    
auto.OR.tree(result.Ser_0.1pTMT.STAR.OR,"0.1% TMT") 
auto.OR.tree(result.Ser_1pTMT.STAR.OR,"1% TMT")
auto.OR.tree(result.Ser_100pTMT.STAR.OR, "100% TMT")
dev.off()

#####################################################


MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
MyGray2 <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
MyGray3 <- rgb(t(col2rgb("black")), alpha=10, maxColorValue=255)
MyTransparent <- rgb(t(col2rgb("black")), alpha=0, maxColorValue=255)
MyGreen <- rgb(t(col2rgb("green")), alpha=150, maxColorValue=255)
MyCyan <- rgb(t(col2rgb("cyan")), alpha=200, maxColorValue=255)
MyRed <- rgb(t(col2rgb("red")), alpha=150, maxColorValue=255)
MyBlue <- rgb(t(col2rgb("blue")), alpha=150, maxColorValue=255)
MyPurple <- rgb(t(col2rgb("purple")), alpha=150, maxColorValue=255)

pdf("05.27.19_ACT.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
result.Ser_0.01pACT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_0.01pACT.STAR,add=1,cpm=0, ttl="0.01p.AcetoPh")
result.Ser_0.1pACT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_0.1pACT.STAR,add=1,cpm=0, ttl="0.1p.AcetoPh")
result.Ser_1pACT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_1pACT.STAR,add=1,cpm=0, ttl="1p.AcetoPh")
result.Ser_100pACT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_100pACT.STAR,add=1,cpm=0, ttl="100p.AcetoPh")
dev.off()

pdf("05.27.19_TMT.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
result.Ser_0.01pTMT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_0.01pTMT.STAR,add=1,cpm=0, ttl="0.01p.TMT")
result.Ser_0.1pTMT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_0.1pTMT.STAR,add=1,cpm=0, ttl="0.1p.TMT")
result.Ser_1pTMT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_1pTMT.STAR,add=1,cpm=0, ttl="1p.TMT")
result.Ser_100pTMT.STAR.OR <- filter.OR.genename(nm.tresult.Ser_100pTMT.STAR,add=1,cpm=0, ttl="100p.TMT")
dev.off()

pdf("07.04.19_ACT_volcano.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
result.Ser_0.01pACT.STAR.OR <- filter.OR.label(nm.result.Ser_0.01pACT.STAR,add=1,cpm=0, ttl="0.01p.AcetoPh")
result.Ser_0.1pACT.STAR.OR <- filter.OR.label(nm.tresult.Ser_0.1pACT.STAR,add=1,cpm=0, ttl="0.1p.AcetoPh")
result.Ser_1pACT.STAR.OR <- filter.OR.label(nm.tresult.Ser_1pACT.STAR,add=1,cpm=0, ttl="1p.AcetoPh")
result.Ser_100pACT.STAR.OR <- filter.OR.label(nm.tresult.Ser_100pACT.STAR,add=1,cpm=0, ttl="100p.AcetoPh")
dev.off()

pdf("07.04.19_TMT_volcano.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
result.Ser_0.01pTMT.STAR.OR <- filter.OR.label(nm.result.Ser_0.01pTMT.STAR,add=1,cpm=0, ttl="0.01p.TMT")
result.Ser_0.1pTMT.STAR.OR <- filter.OR.label(nm.tresult.Ser_0.1pTMT.STAR,add=1,cpm=0, ttl="0.1p.TMT")
result.Ser_1pTMT.STAR.OR <- filter.OR.label(nm.tresult.Ser_1pTMT.STAR,add=1,cpm=0, ttl="1p.TMT")
result.Ser_100pTMT.STAR.OR <- filter.OR.label(nm.tresult.Ser_100pTMT.STAR,add=1,cpm=0, ttl="100p.TMT")
dev.off()

pdf("07.04.19_0.1pTMT.volcano.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
#yaxis height = 70
result.Ser_0.1pTMT.STAR.OR <- filter.OR.label(nm.tresult.Ser_0.1pTMT.STAR,add=1,cpm=0, ttl="0.1p.TMT")
dev.off()          

pdf("06.03.19_ACT_IvanOnlyB.label.pdf")
par(mar=c(6.8,7.3,0.5,0), mgp = c(4.3,1.5,0))
#yaxis height = 40
result.Ser_1pACT.STAR.OR <- filter.OR.label(nm.tresult.Ser_1pACT.STAR,add=1,cpm=0, ttl="1p.AcetoPh")
result.Ser_100pACT.STAR.OR <- filter.OR.label(nm.tresult.Ser_100pACT.STAR,add=1,cpm=0, ttl="100p.AcetoPh")
dev.off()  


filter.OR.label <- function(result,cpm=0,add=0.1,ttl="") { #filter ORs by CPM, correct FDR 
  p1=0.05;p2=0.001;p3=0.001
  result.OR <- subset(result,grepl('Olfr',Gene_name)&(!grepl('-ps',Gene_name)))#|grepl('Taar',Gene_name)|grepl('Vmn',Gene_name))
  result.nonOR <- subset(result,!grepl('Olfr',Gene_name))
  A.cut <- cpm*sum(result$baseMeanA)/1000000; B.cut <- cpm*sum(result$baseMeanB)/1000000
  cut <- cpm*sum(result$baseMean)/1000000
  #  A.cut <- B.cut <- cpm
  sub.result.OR <- subset(result.OR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  sub.result.nonOR <- subset(result.nonOR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  EGR1 <- subset(result.nonOR, Gene_name=="Egr1")
  CFOS <- subset(result.nonOR, Gene_name=="Fos")
  M72 <- subset(result.OR,Gene_name=="Olfr160")
  plot(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),col=MyGray,pch=16)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),bg=MyRed,col=MyGray2,pch=21)
  points(log10(EGR1$baseMeanA+add),log10(EGR1$baseMeanB+add),bg="cyan",pch=21);points(log10(CFOS$baseMeanA+add),log10(CFOS$baseMeanB+add),bg="green",pch=21)
  abline(0,1,lwd=3,lty=2,col="gray")
  legend("bottomright",c("ORs","c-Fos","Egr1","Others"),pt.bg=c(MyRed,"green","cyan",MyTransparent),col=c("black","black","black",MyGray),pch=c(21,21,21,16),bty="n",cex=1.2)
  sub.result.OR$FDR <- p.adjust(sub.result.OR$PValue,method="fdr")
  p1sub <- subset(sub.result.OR,(FDR>p2&FDR<p1)&logFC>0); p2sub <- subset(sub.result.OR,(FDR>p3&FDR<p2)&logFC>0); p3sub <- subset(sub.result.OR,FDR<p3&logFC>0)
  title(ttl)
  MyGray <- "gray"
  plot(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA, cex.axis=1.5)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=1.2)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  #mgp – A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(log10(p1sub$baseMeanA+add),log10(p1sub$baseMeanB+add),col=MyGray2,bg=MyBlue,pch=21,cex=1.2)
  points(log10(p2sub$baseMeanA+add),log10(p2sub$baseMeanB+add),col=MyGray2,bg=MyPurple,pch=21,cex=1.2)
  points(log10(p3sub$baseMeanA+add),log10(p3sub$baseMeanB+add),col=MyGray2,bg=MyRed,pch=21,cex=1.2)
  abline(0,1,lwd=3,lty=2)
  title(ttl)
  lgnd <- c(bquote(italic(p) < .(p3)*' '*(ORs == .(dim(subset(sub.result.OR,FDR<p3&logFC>0))[1]))),
            bquote(italic(p) < .(p1)*' '*(ORs == .(dim(subset(sub.result.OR,FDR<p1&logFC>0))[1]))))
  legend("bottomright", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
  
  #plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA, xlim=c(-4,6.5), ylim=c(0,40), cex.axis=2.5, las=1) #y-axis up to 40
  plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA, xlim=c(-4,6.5), ylim=c(0,70), cex.axis=2.5, las=1) #y-axis up to 70
  #plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(sub.result.OR$logFC,-log10(sub.result.OR$FDR),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=2)
  title(xlab=expression('Log'[2]*'(fold change)'), ylab=expression('-Log'[10]*'(p-value)') ,cex.lab=2.5,mgp=c(4, 1, 0))
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(p1sub$logFC,-log10(p1sub$FDR),col=MyGray2,bg=MyBlue,pch=21,cex=2)
  points(p2sub$logFC,-log10(p2sub$FDR),col=MyGray2,bg=MyPurple,pch=21,cex=2)
  points(p3sub$logFC,-log10(p3sub$FDR),col=MyGray2,bg=MyRed,pch=21,cex=2)
  
  # Olfr923 <- subset(sub.result.OR,Gene_name=="Olfr923")
  # Olfr983 <- subset(sub.result.OR,Gene_name=="Olfr983")
  # Olfr898 <- subset(sub.result.OR,Gene_name=="Olfr898")
  # Olfr888 <- subset(sub.result.OR,Gene_name=="Olfr888")
  # Olfr30 <- subset(sub.result.OR,Gene_name=="Olfr30")
  # Olfr376 <- subset(sub.result.OR,Gene_name=="Olfr376")
  
  # Olfr943<- subset(sub.result.OR,Gene_name=="Olfr943")
  # Olfr933<- subset(sub.result.OR,Gene_name=="Olfr933")
  # Olfr870<- subset(sub.result.OR,Gene_name=="Olfr870")
  # Olfr750<- subset(sub.result.OR,Gene_name=="Olfr750")
  # Olfr745<- subset(sub.result.OR,Gene_name=="Olfr745")
  # Olfr736<- subset(sub.result.OR,Gene_name=="Olfr736")
  # Olfr556<- subset(sub.result.OR,Gene_name=="Olfr556")
  # Olfr490<- subset(sub.result.OR,Gene_name=="Olfr490")
  # Olfr476<- subset(sub.result.OR,Gene_name=="Olfr476")
  # Olfr47<- subset(sub.result.OR,Gene_name=="Olfr47")
  # Olfr452<- subset(sub.result.OR,Gene_name=="Olfr452")
  # Olfr378<- subset(sub.result.OR,Gene_name=="Olfr378")
  # Olfr366<- subset(sub.result.OR,Gene_name=="Olfr366")
  
  # Olfr229<- subset(sub.result.OR,Gene_name=="Olfr229")
  # Olfr1501<- subset(sub.result.OR,Gene_name=="Olfr1501")
  # Olfr1484<- subset(sub.result.OR,Gene_name=="Olfr1484")
  # Olfr1459<- subset(sub.result.OR,Gene_name=="Olfr1459")
  # Olfr142<- subset(sub.result.OR,Gene_name=="Olfr142")
  # Olfr1350<- subset(sub.result.OR,Gene_name=="Olfr1350")
  # Olfr1269<- subset(sub.result.OR,Gene_name=="Olfr1269")
  # Olfr1231<- subset(sub.result.OR,Gene_name=="Olfr1231")
  # Olfr114<- subset(sub.result.OR,Gene_name=="Olfr114")
  # Olfr1135<- subset(sub.result.OR,Gene_name=="Olfr1135")
  # Olfr109<- subset(sub.result.OR,Gene_name=="Olfr109")
  # Olfr1085<- subset(sub.result.OR,Gene_name=="Olfr1085")
  # Olfr1014<- subset(sub.result.OR,Gene_name=="Olfr1014")
  
  # Olfr1395 <- subset(sub.result.OR,Gene_name=="Olfr1395")
  # Olfr1297 <- subset(sub.result.OR,Gene_name=="Olfr1297")
  # Olfr165 <- subset(sub.result.OR,Gene_name=="Olfr165")
  # Olfr135 <- subset(sub.result.OR,Gene_name=="Olfr135")
  # Olfr837 <- subset(sub.result.OR,Gene_name=="Olfr837")
  # Olfr531 <- subset(sub.result.OR,Gene_name=="Olfr531")
  # Olfr376 <- subset(sub.result.OR,Gene_name=="Olfr376")
  
  # Olfr746 <- subset(sub.result.OR,Gene_name=="Olfr736")
  # Olfr556 <- subset(sub.result.OR,Gene_name=="Olfr556")
  # Olfr476 <- subset(sub.result.OR,Gene_name=="Olfr476")
  # Olfr1501 <- subset(sub.result.OR,Gene_name=="Olfr1501")
  # Olfr922 <- subset(sub.result.OR,Gene_name=="Olfr922")
  # Olfr229 <- subset(sub.result.OR,Gene_name=="Olfr229")
  # Olfr47 <- subset(sub.result.OR,Gene_name=="Olfr47")
  # Olfr109 <- subset(sub.result.OR,Gene_name=="Olfr109")
  
  # points(Olfr923$logFC,-log10(Olfr923$FDR),col="gold",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr983$logFC,-log10(Olfr983$FDR),col="darkorange",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr898$logFC,-log10(Olfr898$FDR),col="aquamarine",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr888$logFC,-log10(Olfr888$FDR),col="olivedrab1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr30$logFC,-log10(Olfr30$FDR),col="mediumpurple1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr376$logFC,-log10(Olfr376$FDR),col="mediumorchid1",bg=MyGray3,pch=21,cex=3, lwd=4)
  
  # points(Olfr943$logFC,-log10(Olfr943$FDR),col="gold",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr933$logFC,-log10(Olfr933$FDR),col="darkorange",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr870$logFC,-log10(Olfr870$FDR),col="aquamarine",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr750$logFC,-log10(Olfr750$FDR),col="olivedrab1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr745$logFC,-log10(Olfr745$FDR),col="mediumpurple1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr736$logFC,-log10(Olfr736$FDR),col="mediumorchid1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr556$logFC,-log10(Olfr556$FDR),col="bisque",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr490$logFC,-log10(Olfr490$FDR),col="chartreuse",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr476$logFC,-log10(Olfr476$FDR),col="brown",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr47$logFC,-log10(Olfr47$FDR),col="yellow",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr452$logFC,-log10(Olfr452$FDR),col="cyan3",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr378$logFC,-log10(Olfr378$FDR),col="deepskyblue1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr366$logFC,-log10(Olfr366$FDR),col="pink",bg=MyGray3,pch=21,cex=3, lwd=4)
  
  # points(Olfr229$logFC,-log10(Olfr229$FDR),col="gold",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1501$logFC,-log10(Olfr1501$FDR),col="darkorange",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1484$logFC,-log10(Olfr1484$FDR),col="aquamarine",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1459$logFC,-log10(Olfr1459$FDR),col="olivedrab1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr142$logFC,-log10(Olfr142$FDR),col="mediumpurple1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1350$logFC,-log10(Olfr1350$FDR),col="mediumorchid1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1269$logFC,-log10(Olfr1269$FDR),col="bisque",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1231$logFC,-log10(Olfr1231$FDR),col="chartreuse",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr114$logFC,-log10(Olfr114$FDR),col="brown",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1135$logFC,-log10(Olfr1135$FDR),col="yellow",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr109$logFC,-log10(Olfr109$FDR),col="cyan3",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1085$logFC,-log10(Olfr1085$FDR),col="deepskyblue1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1014$logFC,-log10(Olfr1014$FDR),col="pink",bg=MyGray3,pch=21,cex=3, lwd=4)
  
  # points(Olfr1395$logFC,-log10(Olfr1395$FDR),col="gold",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1297$logFC,-log10(Olfr1297$FDR),col="darkorange",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr165$logFC,-log10(Olfr165$FDR),col="aquamarine",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr135$logFC,-log10(Olfr135$FDR),col="olivedrab1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr837$logFC,-log10(Olfr837$FDR),col="mediumpurple1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr531$logFC,-log10(Olfr531$FDR),col="mediumseagreen",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr376$logFC,-log10(Olfr376$FDR),col="mediumorchid1",bg=MyGray3,pch=21,cex=3, lwd=4)
  
  # points(Olfr746$logFC,-log10(Olfr746$FDR),col="gold",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr556$logFC,-log10(Olfr556$FDR),col="aquamarine",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr476$logFC,-log10(Olfr476$FDR),col="mediumpurple1",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr1501$logFC,-log10(Olfr1501$FDR),col="bisque",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr922$logFC,-log10(Olfr922$FDR),col="chartreuse",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr229$logFC,-log10(Olfr229$FDR),col="brown",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr47$logFC,-log10(Olfr47$FDR),col="cyan3",bg=MyGray3,pch=21,cex=3, lwd=4)
  # points(Olfr109$logFC,-log10(Olfr109$FDR),col="deepskyblue1",bg=MyGray3,pch=21,cex=3, lwd=4)
  
  legend("topleft", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=2.3,bty="n")
  
  #plot.new();legend("topleft", c(as.expression(lgnd), "Olfr943","Olfr933","Olfr870","Olfr750","Olfr745","Olfr736","Olfr556","Olfr490","Olfr476","Olfr47","Olfr452","Olfr378","Olfr366"),pch=21,col=c(MyGray2,MyGray2,"gold","darkorange","aquamarine","olivedrab1","mediumpurple1","mediumorchid1","bisque","chartreuse","brown","yellow","cyan3","deepskyblue1","pink"),pt.bg=c(MyRed,MyBlue,"White","White","White","White","White","White","White","White","White","White","White","White","White"),bty="n",pt.lwd=2,pt.cex=1,cex=1) 
  
  # plot.new();legend("topleft", c(as.expression(lgnd), "Olfr229","Olfr1501","Olfr1484","Olfr1459","Olfr142","Olfr1350","Olfr1269","Olfr1231","Olfr114","Olfr1135","Olfr109","Olfr1085","Olfr1014"),pch=21,col=c(MyGray2,MyGray2,"gold","darkorange","aquamarine","olivedrab1","mediumpurple1","mediumorchid1","bisque","chartreuse","brown","yellow","cyan3","deepskyblue1","pink"),pt.bg=c(MyRed,MyBlue,"White","White","White","White","White","White","White","White","White","White","White","White","White"),bty="n",pt.lwd=2,pt.cex=1,cex=1) 
  
  #legend("topleft", c(as.expression(lgnd), "Olfr923","Olfr983","Olfr30","Olfr376"),col=c(MyGray2,MyGray2,"gold","aquamarine","mediumpurple1","bisque"),pt.bg=c(MyRed,MyBlue,"White","White","White","White"),pch=21,cex=2.2,lwd=c(1,1,4,4,4,4), bty="n")
  #legend("topleft", c(as.expression(lgnd), "Olfr1395","Olfr1297","Olfr165","Olfr376"),col=c(MyGray2,MyGray2,"gold","aquamarine","mediumpurple1","bisque"),pt.bg=c(MyRed,MyBlue,"White","White","White","White"),pch=21,cex=1.5,lwd=c(1,1,4,4,4,4), bty="n")
  #title(ttl)
  return (sub.result.OR)
}

##############

filter.OR.genename <- function(result,cpm=0,add=0.1,ttl="") { #filter ORs by CPM, correct FDR 
  p1=0.05;p2=0.001;p3=0.001
  result.OR <- subset(result,grepl('Olfr',Gene_name)&(!grepl('-ps',Gene_name)))#|grepl('Taar',Gene_name)|grepl('Vmn',Gene_name))
  result.nonOR <- subset(result,!grepl('Olfr',Gene_name))
  A.cut <- cpm*sum(result$baseMeanA)/1000000; B.cut <- cpm*sum(result$baseMeanB)/1000000
  cut <- cpm*sum(result$baseMean)/1000000
  #  A.cut <- B.cut <- cpm
  sub.result.OR <- subset(result.OR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  sub.result.nonOR <- subset(result.nonOR, baseMeanA>=A.cut&baseMeanB>=B.cut)
  EGR1 <- subset(result.nonOR, Gene_name=="Egr1")
  CFOS <- subset(result.nonOR, Gene_name=="Fos")
  M72 <- subset(result.OR,Gene_name=="Olfr160")
  plot(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(log10(sub.result.nonOR$baseMeanA+add),log10(sub.result.nonOR$baseMeanB+add),col=MyGray,pch=16)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),bg=MyRed,col=MyGray2,pch=21)
  points(log10(EGR1$baseMeanA+add),log10(EGR1$baseMeanB+add),bg="cyan",pch=21);points(log10(CFOS$baseMeanA+add),log10(CFOS$baseMeanB+add),bg="green",pch=21)
  abline(0,1,lwd=3,lty=2,col="gray")
  legend("bottomright",c("ORs","c-Fos","Egr1","Others"),pt.bg=c(MyRed,"green","cyan",MyTransparent),col=c("black","black","black",MyGray),pch=c(21,21,21,16),bty="n",cex=1.2)
  sub.result.OR$FDR <- p.adjust(sub.result.OR$PValue,method="fdr")
  p1sub <- subset(sub.result.OR,(FDR>p2&FDR<p1)&logFC>0); p2sub <- subset(sub.result.OR,(FDR>p3&FDR<p2)&logFC>0); p3sub <- subset(sub.result.OR,FDR<p3&logFC>0)
  title(ttl)
  
  MyGray <- "gray"
  plot(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),type="n",frame.plot=F,xlab=NA,ylab=NA, cex.axis=1.5)
  points(log10(sub.result.OR$baseMeanA+add),log10(sub.result.OR$baseMeanB+add),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=1.2)
  title(xlab=expression('Log'[10]*'(mean read counts, no-odor control)'), ylab=expression('Log'[10]*'(mean read counts, stimulated)') ,cex.lab=1.5)
  #mgp – A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(log10(p1sub$baseMeanA+add),log10(p1sub$baseMeanB+add),col=MyGray2,bg=MyBlue,pch=21,cex=1.2)
  points(log10(p2sub$baseMeanA+add),log10(p2sub$baseMeanB+add),col=MyGray2,bg=MyPurple,pch=21,cex=1.2)
  points(log10(p3sub$baseMeanA+add),log10(p3sub$baseMeanB+add),col=MyGray2,bg=MyRed,pch=21,cex=1.2)
  abline(0,1,lwd=3,lty=2)
  title(ttl)
  lgnd <- c(bquote(italic(p) < .(p3)*' '*(ORs == .(dim(subset(sub.result.OR,FDR<p3&logFC>0))[1]))),
            bquote(italic(p) < .(p1)*' '*(ORs == .(dim(subset(sub.result.OR,FDR<p1&logFC>0))[1]))))
  legend("bottomright", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=1.2,bty="n")
  MyGray <- rgb(t(col2rgb("black")), alpha=60, maxColorValue=255)
  #plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA, xlim=c(-5,5), ylim=c(0,70), cex.axis=2.5, las=1) #ACT y-axis up to 40
  plot(sub.result.OR$logFC,-log10(sub.result.OR$FDR),type="n",frame.plot=F,xlab=NA,ylab=NA)
  points(sub.result.OR$logFC,-log10(sub.result.OR$FDR),
         bg=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),
         col=ifelse(sub.result.OR$FDR<p1&sub.result.OR$logFC>0, MyTransparent, MyGray),pch=21,cex=2)
  title(xlab=expression('Log'[2]*'(fold change)'), ylab=expression('-Log'[10]*'(p-value)') ,cex.lab=2.5,mgp=c(4, 1, 0))
  #title(xlab="Mean read counts, no-odor control, log10", ylab="Mean read counts, stimulated, log10",cex.lab=1.5)
  points(p1sub$logFC,-log10(p1sub$FDR),col=MyGray2,bg=MyBlue,pch=21,cex=2)
  points(p2sub$logFC,-log10(p2sub$FDR),col=MyGray2,bg=MyPurple,pch=21,cex=2)
  points(p3sub$logFC,-log10(p3sub$FDR),col=MyGray2,bg=MyRed,pch=21,cex=2)
  
  legend("topleft", as.expression(lgnd),,col=c(MyGray2,MyGray2),pt.bg=c(MyRed,MyBlue),pch=21,cex=2.3,bty="n")
  
  return (sub.result.OR)
}          

ACT.0.01p.Taar <- subset(nm.result.Ser_0.01pACT.STAR,grepl('Taar',Gene_name))
TMT.0.01p.Taar <- subset(nm.result.Ser_0.01pTMT.STAR,grepl('Taar',Gene_name))
ACT.0.1p.Taar <- subset(nm.tresult.Ser_0.1pACT.STAR,grepl('Taar',Gene_name))
TMT.0.1p.Taar <- subset(nm.tresult.Ser_0.1pTMT.STAR,grepl('Taar',Gene_name))
ACT.100p.Taar <- subset(nm.tresult.Ser_100pACT.STAR,grepl('Taar',Gene_name))
TMT.100p.Taar <- subset(nm.tresult.Ser_100pTMT.STAR,grepl('Taar',Gene_name))
ACT.1p.Taar <- subset(nm.tresult.Ser_1pACT.STAR,grepl('Taar',Gene_name))
TMT.1p.Taar <- subset(nm.tresult.Ser_1pTMT.STAR,grepl('Taar',Gene_name))
write.csv(ACT.0.01p.Taar, file = "ACT.0.01p.Taar.csv")
write.csv(TMT.0.01p.Taar, file = "TMT.0.01p.Taar.csv")
write.csv(ACT.0.1p.Taar, file = "ACT.0.1p.Taar.csv")
write.csv(TMT.0.1p.Taar, file = "TMT.0.1p.Taar.csv")
write.csv(ACT.100p.Taar, file = "ACT.100p.Taar.csv")
write.csv(TMT.100p.Taar, file = "TMT.100p.Taar.csv")
write.csv(ACT.1p.Taar, file = "ACT.1p.Taar.csv")
write.csv(TMT.1p.Taar, file = "TMT.1p.Taar.csv")

#Differential Expression function
DiffExp <- function (targets, countsTable) {
  Treat <- factor(targets$Treatment);Subject <- factor(targets$Subject);design <- model.matrix(~Subject+Treat)
  cds <- newCountDataSet(countsTable,Treat);cds <- estimateSizeFactors(cds);cds <- estimateDispersions(cds);d <- nbinomTest(cds,"0","1")
  e.litter <- DGEList(counts=countsTable)
  e.litter <- estimateGLMCommonDisp(e.litter,design)
  e.litter <- estimateGLMTrendedDisp(e.litter,design)
  e.litter <- estimateGLMTagwiseDisp(e.litter,design)
  fit <- glmFit(e.litter, design);lrt <- glmLRT(fit);diff <- topTags(lrt,n=dim(lrt)[1])$table
  result <- merge(merge(diff,countsTable,by=0,sort=F),d, by.x="Row.names", by.y="id",sort=F)
  colnames(result)[1] <- "id"
  return(result)
}

targets.ST3 <- data.frame(FileName=c("no1","no2","no3","ACT1","ACT2","ACT3"),Subject=c("A","B","C","A","B","C"),Treatment=c(0,0,0,1,1,1))

targets.ST4 <- data.frame(FileName=c("No1","No5","No6","No7","TMT1","TMT5","TMT6","TMT7"),Subject=c("A","B","C","D","A","B","C","D"),Treatment=c(0,0,0,0,1,1,1,1))

#DREAM 5% ACT Data Ivan
Ivan1R <- read.table("1R_X2l5xLrWCivJ_L1_R1_001.genes.results", header = T) #Cont
Ivan2R <- read.table("2R_4T0AGN3UcppS_L1_R1_001.genes.results", header = T) #Cont
Ivan3R <- read.table("3R_ajiAh4sp6fIG_L1_R1_001.genes.results", header = T) #Cont
Ivan10R <- read.table("10R_re2F5zjmGZ7e_L2_R1_001.genes.results", header = T) #IvanAceto
Ivan11R <- read.table("11R_Qkx50qMje0Ic_L2_R1_001.genes.results", header = T) #IvanAceto
Ivan12R <- read.table("12R_YVRhY3k5ZXpZ_L2_R1_001.genes.results", header = T) #IvanAceto
targets.ST3 <- data.frame(FileName=c("no1","no2","no3","ACT1","ACT2","ACT3"),Subject=c("A","B","C","A","B","C"),Treatment=c(0,0,0,1,1,1))
countsTable.Ivan_5pACT.STAR <- round(data.frame(no1=Ivan1R$expected_count,
                                                no2=Ivan2R$expected_count,
                                                no3=Ivan3R$expected_count,
                                                st1=Ivan10R$expected_count,
                                                st2=Ivan11R$expected_count,
                                                st3=Ivan12R$expected_count,
                                                row.names=Ivan12R$gene_id))
result.Ivan_5pACT.STAR  <- DiffExp(targets.ST3,countsTable.Ivan_5pACT.STAR)
tresult.Ivan_5pACT.STAR <- as_tibble(result.Ivan_5pACT.STAR) %>% rename("gene_id" = "id")
nm.tresult.Ivan_5pACT.STAR <- left_join(tresult.Ivan_5pACT.STAR, ens_names, by = "gene_id") %>% select(Gene_name, logFC:padj) %>% mutate(isOlfr = ifelse(str_detect(Gene_name, "Olfr"), "isOlfr", "notOlfr")) %>% mutate(isSignif = ifelse(FDR < 0.05, "Pval<0.05", "Notsig"))
write.csv(nm.tresult.Ivan_5pACT.STAR, file = "All_result.Ivan_5pACT.STAR.csv")
pdf("result.Ivan_5pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ivan_5pACT.STAR.OR <- filter.OR.label(nm.tresult.Ivan_5pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ivan_5pACT.STAR.OR, file = "OR_result.Ivan_5pACT.STAR.csv")


#Yue 1% ACT Data 
YJ3.STAR <- read.table("YJ3_AGGCAGAA_L008_R1_001.genes.results", header = T) #Cont
YJA1.STAR <- read.table("YJA1_TAAGGCGA_L008_R1.genes.results", header = T) #Cont
YJA4.STAR <- read.table("YJA4_TCCTGAGC_L008_R1.genes.results", header = T) #Cont
YJ4.STAR <- read.table("YJ4_TCCTGAGC_L008_R1_001.genes.results", header = T) #IvanAceto
YJA3.STAR <- read.table("YJA3_AGGCAGAA_L008_R1.genes.results", header = T) #IvanAceto
YJA6.STAR <- read.table("YJA6_TAGGCATG_L008_R1.genes.results", header = T) #IvanAceto
targets.ST3 <- data.frame(FileName=c("no1","no2","no3","ACT1","ACT2","ACT3"),Subject=c("A","B","C","A","B","C"),Treatment=c(0,0,0,1,1,1))
countsTable.Yue_1pACT.STAR <- round(data.frame(no1=YJ3.STAR$expected_count,
                                               no2=YJA1.STAR$expected_count,
                                               no3=YJA4.STAR$expected_count,
                                               st1=YJ4.STAR$expected_count,
                                               st2=YJA3.STAR$expected_count,
                                               st3=YJA6.STAR$expected_count,
                                               row.names=YJA6.STAR$gene_id))
result.Yue_1pACT.STAR  <- DiffExp(targets.ST3,countsTable.Yue_1pACT.STAR)
tresult.Yue_1pACT.STAR <- as_tibble(result.Yue_1pACT.STAR) %>% rename("gene_id" = "id")
nm.tresult.Yue_1pACT.STAR <- left_join(tresult.Yue_1pACT.STAR, ens_names, by = "gene_id") %>% select(Gene_name, logFC:padj) %>% mutate(isOlfr = ifelse(str_detect(Gene_name, "Olfr"), "isOlfr", "notOlfr")) %>% mutate(isSignif = ifelse(FDR < 0.05, "Pval<0.05", "Notsig"))
write.csv(nm.tresult.Yue_1pACT.STAR, file = "All_result.Yue_1pACT.STAR.csv")
pdf("result.Yue_1pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Yue_1pACT.STAR.OR <- filter.OR.label(nm.tresult.Yue_1pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Yue_1pACT.STAR.OR, file = "OR_result.Yue_1pACT.STAR.csv")

#Serene Dor/Ventral Data
YJK4.STAR <- read.table("YJK4_TTAGGC.genes.results", header=T) #Dor
YJK5.STAR <- read.table("YJK5_GATCAG.genes.results", header=T) #Ven
YJK6.STAR <- read.table("YJK6_GTTTCG.genes.results", header=T) #Dor
YJK7.STAR <- read.table("YJK7_CGTACG.genes.results", header=T) #Ven
YJK8.STAR <- read.table("YJK8_GAGTGG.genes.results", header=T) #Dor
YJK9.STAR <- read.table("YJK9_ATTCCT.genes.results", header=T) #Ven

#Serene 0.01p ACT Data
SerE1.STAR <- read.table("SerE1_L001_R1_001.genes.results", header=T) #no
SerE4.STAR <- read.table("SerE4_L001_R1_001.genes.results", header=T) #no
SerE7.STAR <- read.table("SerE7_L001_R1_001.genes.results", header=T) #no
SerE2.STAR <- read.table("SerE2_L001_R1_001.genes.results", header=T)#ACT
SerE5.STAR <- read.table("SerE5_L001_R1_001.genes.results", header=T)#ACT
SerE8.STAR <- read.table("SerE8_L001_R1_001.genes.results", header=T)#ACT
#Serene 0.1p ACT Data
SerE3.STAR <- read.table("SerE3_L001_R1_001.genes.results", header=T)#ACT
SerE6.STAR <- read.table("SerE6_L001_R1_001.genes.results", header=T)#ACT
SerE9.STAR <- read.table("SerE9_L001_R1_001.genes.results", header=T)#ACT

#Serene 0.01p TMT Data
SerG1.STAR <- read.table("SerG1_L001_R1_001.genes.results", header=T) #no
SerG5.STAR <- read.table("SerG5_L001_R1_001.genes.results", header=T) #no
SerG9.STAR <- read.table("SerG9_L001_R1_001.genes.results", header=T) #no
SerG2.STAR <- read.table("SerG2_L001_R1_001.genes.results", header=T)#TMT
SerG6.STAR <- read.table("SerG6_L001_R1_001.genes.results", header=T)#TMT
SerG10.STAR <- read.table("SerG10_L001_R1_001.genes.results", header=T)#TMT
#Serene 0.1p TMT Data
SerG3.STAR <- read.table("SerG3_L001_R1_001.genes.results", header=T)#TMT
SerG7.STAR <- read.table("SerG7_L001_R1_001.genes.results", header=T)#TMT
SerG11.STAR <- read.table("SerG11_L001_R1_001.genes.results", header=T)#TMT

countsTable.Yue_DorVen.STAR <- round(data.frame(Dorsal1=YJK4.STAR$expected_count,
                                                Dorsal2=YJK6.STAR$expected_count,
                                                Dorsal3=YJK8.STAR$expected_count,
                                                Ventral1=YJK5.STAR$expected_count,
                                                Ventral2=YJK7.STAR$expected_count,
                                                Ventral3=YJK9.STAR$expected_count,
                                                row.names=YJK9.STAR$gene_id))

targets.DV3 <- data.frame(FileName=c("Dorsal1","Dorsal2","Dorsal3"," Ventral1"," Ventral2"," Ventral3"),Subject=c("A","B","C","A","B","C"),Treatment=c(0,0,0,1,1,1))

countsTable.Ser_0.1pACT.STAR <- round(data.frame(no1=SerE1.STAR$expected_count,
                                                 no2=SerE4.STAR$expected_count,
                                                 no3=SerE7.STAR$expected_count,
                                                 st1=SerE3.STAR$expected_count,
                                                 st2=SerE6.STAR$expected_count,
                                                 st3=SerE9.STAR$expected_count,
                                                 row.names=SerE9.STAR$gene_id))

countsTable.Ser_0.1pTMT.STAR <- round(data.frame(no1=SerG1.STAR$expected_count,
                                                 no2=SerG5.STAR$expected_count,
                                                 no3=SerG9.STAR$expected_count,
                                                 st1=SerG3.STAR$expected_count,
                                                 st2=SerG7.STAR$expected_count,
                                                 st3=SerG11.STAR$expected_count,
                                                 row.names=SerG11.STAR$gene_id))

result.Yue_DorVen.STAR  <- DiffExp(targets.DV3,countsTable.Yue_DorVen.STAR)
tresult.Yue_DorVen.STAR <- as_tibble(result.Yue_DorVen.STAR) %>% rename("gene_id" = "id")
nm.tresult.Yue_DorVen.STAR <- left_join(tresult.Yue_DorVen.STAR, ens_names, by = "gene_id") %>% select(Gene_name, logFC:padj) %>% mutate(isOlfr = ifelse(str_detect(Gene_name, "Olfr"), "isOlfr", "notOlfr")) %>% mutate(isSignif = ifelse(FDR < 0.05, "Pval<0.05", "Notsig"))
write.csv(nm.tresult.Yue_DorVen.STAR, file = "All_result.Yue_DorVen.STAR.csv")
pdf("result.Yue_DorVen.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Yue_DorVen.STAR.OR <- filter.OR.genename(nm.tresult.Yue_DorVen.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Yue_DorVen.STAR.OR, file = "OR_result.Yue_DorVen.STAR.csv")

result.Ser_0.1pACT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_0.1pACT.STAR)
write.csv(result.Ser_0.1pACT.STAR, file = "All_result.Ser_0.1pACT.STAR.csv")
pdf("result.Ser_0.1pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_0.1pACT.STAR.OR <- filter.OR (result.Ser_0.1pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_0.1pACT.STAR.OR, file = "OR_result.Ser_0.1pACT.STAR.csv")

result.Ser_0.1pTMT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_0.1pTMT.STAR)
write.csv(result.Ser_0.1pTMT.STAR, file = "All_result.Ser_0.1pTMT.STAR.csv")
pdf("result.Ser_0.1pTMT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_0.1pTMT.STAR.OR <- filter.OR (result.Ser_0.1pTMT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_0.1pTMT.STAR.OR, file = "OR_result.Ser_0.1pTMT.STAR.csv")


countsTable.Ser_0.01pACT.STAR <- round(data.frame(no1=SerE1.STAR$expected_count,
                                                  no2=SerE4.STAR$expected_count,
                                                  no3=SerE7.STAR$expected_count,
                                                  st1=SerE2.STAR$expected_count,
                                                  st2=SerE5.STAR$expected_count,
                                                  st3=SerE8.STAR$expected_count,
                                                  row.names=SerE8.STAR$gene_id))

countsTable.Ser_0.01pTMT.STAR <- round(data.frame(no1=SerG1.STAR$expected_count,
                                                  no2=SerG5.STAR$expected_count,
                                                  no3=SerG9.STAR$expected_count,
                                                  st1=SerG2.STAR$expected_count,
                                                  st2=SerG6.STAR$expected_count,
                                                  st3=SerG10.STAR$expected_count,
                                                  row.names=SerG10.STAR$gene_id))

result.Ser_0.01pACT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_0.01pACT.STAR)
write.csv(result.Ser_0.01pACT.STAR, file = "All_result.Ser_0.01pACT.STAR.csv")
pdf("result.Ser_0.01pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_0.01pACT.STAR.OR <- filter.OR (result.Ser_0.01pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_0.01pACT.STAR.OR, file = "OR_result.Ser_0.01pACT.STAR.csv")

result.Ser_0.01pTMT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_0.01pTMT.STAR)
write.csv(result.Ser_0.01pTMT.STAR, file = "All_result.Ser_0.01pTMT.STAR.csv")
pdf("result.Ser_0.01pTMT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_0.01pTMT.STAR.OR <- filter.OR (result.Ser_0.01pTMT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_0.01pTMT.STAR.OR, file = "OR_result.Ser_0.01pTMT.STAR.csv")

#Serene 100p ACT Data#
SerA4.STAR <- read.table("SerA4_TCCTGAGC.genes.results", header=T) #no
SerB4.STAR <- read.table("SerB4_TCCTGAGC.genes.results", header=T) #no
SerB10.STAR <- read.table("SerB10_CGAGGCTG.genes.results", header=T) #no
SerA6.STAR <- read.table("SerA6_TAGGCATG.genes.results", header=T) #ACT
SerB6.STAR <- read.table("SerB6_TAGGCATG.genes.results", header=T) #ACT
SerB12.STAR <- read.table("SerB12_GTAGAGGA.genes.results", header=T) #ACT
#Serene 1p ACT Data
SerA5.STAR <- read.table("SerA5_GGACTCCT.genes.results", header=T) #ACT
SerB5.STAR <- read.table("SerB5_GGACTCCT.genes.results", header=T) #ACT
SerB11.STAR <- read.table("SerB11_AAGAGGCA.genes.results", header=T) #ACT

#Serene 100p TMT Data#
SerA7.STAR <- read.table("SerA7_CTCTCTAC.genes.results", header=T) #no
SerA10.STAR <- read.table("SerA10_CGAGGCTG.genes.results", header=T) #no
SerB1.STAR <- read.table("SerB1_TAAGGCGA.genes.results", header=T) #no
SerA9.STAR <- read.table("SerA9_GCTACGCT.genes.results", header=T) #TMT
SerA12.STAR <- read.table("SerA12_GTAGAGGA.genes.results", header=T) #TMT
SerB3.STAR <- read.table("SerB3_AGGCAGAA.genes.results", header=T) #TMT
#Serene 1p TMT Data
SerA8.STAR <- read.table("SerA8_CAGAGAGG.genes.results", header=T) #TMT
SerA11.STAR <- read.table("SerA11_AAGAGGCA.genes.results", header=T) #TMT
SerB2.STAR <- read.table("SerB2_CGTACTAG.genes.results", header=T) #TMT

countsTable.Ser_100pACT.STAR <- round(data.frame(no1=SerA4.STAR$expected_count,
                                                 no2=SerB4.STAR$expected_count,
                                                 no3=SerB10.STAR$expected_count,
                                                 st1=SerA6.STAR$expected_count,
                                                 st2=SerB6.STAR$expected_count,
                                                 st3=SerB12.STAR$expected_count,
                                                 row.names=SerB12.STAR$gene_id))

countsTable.Ser_100pTMT.STAR <- round(data.frame(no1=SerA7.STAR$expected_count,
                                                 no2=SerA10.STAR$expected_count,
                                                 no3=SerB1.STAR$expected_count,
                                                 st1=SerA9.STAR$expected_count,
                                                 st2=SerA12.STAR$expected_count,
                                                 st3=SerB3.STAR$expected_count,
                                                 row.names=SerB3.STAR$gene_id))

result.Ser_100pACT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_100pACT.STAR)
write.csv(result.Ser_100pACT.STAR, file = "All_result.Ser_100pACT.STAR.csv")
pdf("result.Ser_100pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_100pACT.STAR.OR <- filter.OR (result.Ser_100pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_100pACT.STAR.OR, file = "OR_result.Ser_100pACT.STAR.csv")

result.Ser_100pTMT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_100pTMT.STAR)
write.csv(result.Ser_100pTMT.STAR, file = "All_result.Ser_100pTMT.STAR.csv")
pdf("result.Ser_100pTMT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_100pTMT.STAR.OR <- filter.OR (result.Ser_100pTMT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_100pTMT.STAR.OR, file = "OR_result.Ser_100pTMT.STAR.csv")

#########
countsTable.Ser_1pACT.STAR <- round(data.frame(no1=SerA4.STAR$expected_count,
                                               no2=SerB4.STAR$expected_count,
                                               no3=SerB10.STAR$expected_count,
                                               st1=SerA5.STAR$expected_count,
                                               st2=SerB5.STAR$expected_count,
                                               st3=SerB11.STAR$expected_count,
                                               row.names=SerB11.STAR$gene_id))

countsTable.Ser_1pTMT.STAR <- round(data.frame(no1=SerA7.STAR$expected_count,
                                               no2=SerA10.STAR$expected_count,
                                               no3=SerB1.STAR$expected_count,
                                               st1=SerA8.STAR$expected_count,
                                               st2=SerA11.STAR$expected_count,
                                               st3=SerB2.STAR$expected_count,
                                               row.names=SerB2.STAR$gene_id))

result.Ser_1pACT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_1pACT.STAR)
write.csv(result.Ser_1pACT.STAR, file = "All_result.Ser_1pACT.STAR.csv")
pdf("result.Ser_1pACT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_1pACT.STAR.OR <- filter.OR (result.Ser_1pACT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_1pACT.STAR.OR, file = "OR_result.Ser_1pACT.STAR.csv")

result.Ser_1pTMT.STAR  <- DiffExp(targets.ST3,countsTable.Ser_1pTMT.STAR)
write.csv(result.Ser_1pTMT.STAR, file = "All_result.Ser_1pTMT.STAR.csv")
pdf("result.Ser_1pTMT.STAR.pdf",5.5,5.5)
par(mar=c(5.5,5.5,2,1))
result.Ser_1pTMT.STAR.OR <- filter.OR (result.Ser_1pTMT.STAR,add=1,cpm=0)
dev.off() #Closes the pdf
write.csv(result.Ser_1pTMT.STAR.OR, file = "OR_result.Ser_1pTMT.STAR.csv")


