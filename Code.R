
######################1. GEO dataset merge##########################

library(data.table)
library(tidyverse)
library(limma)
library(sva)
library(dplyr)

###Combine 3 GEO datasets (GSE39582，GSE72970，GSE103479) and combat

getwd()
setwd("/Users/qianruzhu/Desktop/HZMUST/1_GEO")

#GSE39582
probeMatrix=read.table("GSE39582_probeMatrix.txt", sep = "\t", header = T)
rownames(probeMatrix)=probeMatrix$ID_REF
ID=read.table("GSE39582_ID.txt", sep="\t", header = T)
probeMatrix=probeMatrix[ID$ID_REF,]
EXP=merge(ID, probeMatrix, by='ID_REF')
EXP=EXP[,-1]
rownames(EXP)=EXP$Gene
index=duplicated(EXP$Gene)
EXP=EXP[!index,]
rownames(EXP)=EXP$Gene
EXP=EXP[,-1]
write.table(EXP, file="GSE39582_EXP.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

#GSE72970
probeMatrix=read.table("GSE79720_probeMatrix.txt", sep = "\t", header = T)
rownames(probeMatrix)=probeMatrix$ID_REF
ID=read.table("GSE79720_ID.txt", sep="\t", header = T)
probeMatrix=probeMatrix[ID$ID_REF,]
EXP=merge(ID, probeMatrix, by='ID_REF')
EXP=EXP[,-1]
rownames(EXP)=EXP$Gene
index=duplicated(EXP$Gene)
EXP=EXP[!index,]
rownames(EXP)=EXP$Gene
EXP=EXP[,-1]
write.table(EXP, file="GSE79720_EXP.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

#GSE103479
probeMatrix=read.table("GSE103479_probeMatrix.txt", sep = "\t", header = T)
rownames(probeMatrix)=probeMatrix$ID_REF
ID=read.table("GSE103479_ID.txt", sep="\t", header = T)
probeMatrix=probeMatrix[ID$ID_REF,]
EXP=merge(ID, probeMatrix, by='ID_REF')
EXP=EXP[,-1]
rownames(EXP)=EXP$Gene
index=duplicated(EXP$Gene)
EXP=EXP[!index,]
rownames(EXP)=EXP$Gene
EXP=EXP[,-1]
write.table(EXP, file="GSE103479_EXP.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

#acombine and combat
setwd("/Users/qianruzhu/Desktop/HZMUST/2_GEOmerge")
outFile="merge.txt"      

#Get all files ending in ".txt" in the directory
files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#Read all files ending in ".txt" in the directory and save to genelist
for(file in files){
  if(file==outFile){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #read input file
  geneNames=as.vector(rt[,1])        #Extract gene names
  uniqGene=unique(geneNames)       #acquire Unique gene 
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#acquire intersect gene
interGenes=Reduce(intersect, geneList)

#merge dataset
exp1 <- read.table("GSE39582.txt", header=T, row.names=1, sep="\t")      
exp2 <- read.table("GSE79720.txt", header=T, row.names=1, sep="\t")
exp3 <- read.table("GSE103479.txt", header=T, row.names=1, sep="\t")     
exp1 <- exp1[interGenes,]
exp2 <- exp2[interGenes,]
exp3 <- exp3[interGenes,]
rt=cbind.data.frame(exp1,exp2,exp3)
batchType=c(rep(1,585),rep(2,124),rep(3,156))

#data correction
outTab=ComBat(rt, batchType, par.prior=TRUE)
outtab=data.frame(outTab)
write.table(outTab, file="merge_all.txt", sep="\t", quote=F, col.names=F)

#Extract data of 804 patients
clinical <- read.table("clinical_804.txt", sep="\t", header = T)
merge <- as.data.frame(t(outTab))
merge_804 <- merge[clinical$Barcode,]
write.table(merge_804, file="merge_804.txt", sep="\t", quote=F, col.names=T, row.names = T)

#----------------------------------------------------Fig1---------------------------------------------------
# Fig1_1.m6Adiff
#Extract the m6A expression matrix
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig1/m6Adiff")
m6A <- read.table("m6A.txt", sep="\t",header=T)
merge <- read.table("merge_804.txt",sep="\t",row.names=1,header=T)
merge <- t(merge)
merge_m6A=as.matrix(merge[m6A$m6A,])

#m6A diff analysis
library(limma)
library(reshape2)
library(ggpubr)

#Convert to ggplot2 input file
exp=merge_m6A
exp=as.data.frame(t(exp))
exp <- log2(exp+1)
sample=as.matrix(rownames(exp))
clinical <- read.table("clinical_804.txt", sep="\t", header = T)
sampleType <- ifelse(clinical$RFSevent=="0", "no recurrence", "recurrence")
group=cbind.data.frame(sample,sampleType)
exp=cbind(exp, Type=group$sampleType)

data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#draw  boxplot
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c("#4169E1","#B22222"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#output 'pdf'
pdf(file="boxplot_recurrence_log2.pdf", width=7, height=5)
print(p1)
dev.off()

#clinical_stage
clinical <- read.table("clinical_stage.txt", sep="\t", header = T)
merge <- read.table("merge_804.txt",sep="\t",row.names=1,header=T)
merge_stage <- merge[clinical$Sample_geo_accession,]


#extract the m6A expression matrix
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig1")
m6A <- read.table("m6A.txt", sep="\t",header=T)
merge <- t(merge_stage)
merge_m6A=as.data.frame(merge[m6A$m6A,])

#m6A diff analysis
library(limma)
library(reshape2)
library(ggpubr)

#Convert to ggplot2 input file
exp=merge_m6A
exp=as.data.frame(t(exp))
exp <- log2(exp+1)
sample=as.matrix(rownames(exp))
sampleType <- ifelse(clinical$TNM=="I-II", "I-II", "III-IV")
group=cbind.data.frame(sample,sampleType)
exp=cbind(exp, Type=group$sampleType)

data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#boxplot
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c("#4169E1","#B22222"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#pdf
pdf(file="boxplot_stage_log2.pdf", width=7, height=5)
print(p1)
dev.off()

#fig1_2.CNVfreq
# read the m6A expression matrix
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig1/CNVfeq")
rt=read.table("CNV_all.txt", header=T, sep="\t", check.names=F)  
rownames(rt)=rt$Gene
index=duplicated(rt$Gene)
rt=rt[!index,]
rownames(rt)=rt$Gene
rt=rt[,-1]
m6A=read.table("m6A.txt", header=T, sep="\t", check.names=F)  
CNV_m6A=rt[m6A$m6A,]
CNV_m6A=t(CNV_m6A)
write.table(CNV_m6A, "CNV_m6A.txt", sep="\t", row.names = T)

clinical=read.table("clinical_RFS.txt", header=T, sep="\t", check.names=F)  
CNV_m6A=read.table("CNV_m6A.txt", header=T, sep="\t", row.names=1,check.names=F) 
samesample <- intersect(clinical$id,row.names(CNV_m6A))
CNV_m6A=CNV_m6A[samesample,]
rt=CNV_m6A
rt=t(rt)
GAIN=rowSums(rt> 0)       #the number of increase
LOSS=rowSums(rt< 0)       #he number of decrease
GAIN=GAIN/ncol(rt)*100      #the percentage of increase
LOSS=LOSS/ncol(rt)*100      #the percentage of decrease
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#draw
data.max = apply(data, 1, max)
pdf(file="CNVfreq_recurrence.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()

#CNVcir

setwd("/Users/qianruzhu/Desktop/HZMUST/Fig1/CNVcir")
install.packages("RCircos")
library(RCircos)
#read file
cnvMatrix=read.table("CNV_m6A_recurrence_Y.txt", sep="\t",row.names = 1, header = T)
cnvMatrix=as.data.frame(t(cnvMatrix))
geneRef=read.table("geneRef.txt", sep="\t", row.names=1,header = T)
geneRef=geneRef[colnames(cnvMatrix),]
geneRef=geneRef[rownames(cnvMatrix),]
write.table(geneRef, "Rcircos.genelable.txt",sep="\t", row.names = T)

#circle
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#Set circle diagram parameters
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=1
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#output
pdf(file="RCircos.pdf", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#scatter plot
RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t",check.names=F)
RCircos.Scatter.Data=as.data.frame(RCircos.Scatter.Data)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

#add gene name
RCircos.Gene.Label.Data=read.table("Rcircos.genelable.txt", header=T, sep="\t", check.names=F)
name.col <- 5
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()
data(UCSC.HG19.Human.CytoBandIdeogram)

#--------------------------------------------------------Fig2-------------------------------------------------------
#Fig2_1:m6A_survival 
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig2/m6Asurvival")
#library
library(limma)
library(survival)
library(survminer)
expFile="m6aGeneExp.txt"    
cliFile="clinical_804.txt"          

#read expression file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#read clinical file
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     
cli$RFS <- cli$RFS/12
cli <- cli[,-3]

#merge clinical and expression file
rt=cbind(cli, data)

#对基因进行循环，找出预后相关基因
outTab=data.frame()
km=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox analysis
  cox <- coxph(Surv(RFS,RFSevent) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  #km analysis
  data=rt[,c("RFS", "RFSevent", i)]
  colnames(data)=c("RFS", "RFSevent", "gene")
  #cutoff
  res.cut=surv_cutpoint(data, time = "RFS", event = "RFSevent", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(RFS,RFSevent) ~gene, data = res.cat)
  #print(paste0(i, " ", res.cut$cutpoint[1]))
  diff=survdiff(Surv(RFS,RFSevent) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  km=c(km, pValue)
  #pvalue<0.05
  if(pValue<0.05){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    
    #output
    surPlot=ggsurvplot(fit,
                       data=res.cat,
                       pval=pValue,
                       pval.size=6,
                       legend.title=i,
                       legend.labs=c("high","low"),
                       xlab="Time(years)",
                       ylab="Recurrence survival",
                       palette=c("#B22222", "#4169E1"), 
                       break.time.by=1,
                       conf.int=T,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("sur.", i, ".pdf"),onefile = FALSE,
        width = 7,        
        height =5)       
    print(surPlot)
    dev.off()
  }
}

#output
outTab=cbind(outTab, km)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
#forest plot
bioForest=function(coxFile=null,forestFile=null){
  #read input file
  rt=read.table("uniCox_forest.txt", header=T, sep="\t", row.names=1, check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrLow[hrLow<0.001]=0.001
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #figure
  height=(nrow(rt)/15)+5.5
  pdf(file=forestFile, width=7, height=height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #Plot the clinical information on the left side of the forest plot
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  #draw forest plot
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  LOGindex=2 
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "red", "blue")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,LOGindex^a1)
  dev.off()
}
bioForest(coxFile="uniCox_forest.txt", forestFile="forest_m6A_cox.pdf")

# m6A network
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig2/m6Anetwork")
install.packages("igraph")
install.packages("psych")
install.packages("reshape2")
install.packages("RColorBrewer")

library(igraph)
library(psych)
library(reshape2)
library("RColorBrewer")
GeneExpfile <- "m6aGeneExp.txt"    
Genefile <- "gene.txt"              
Coxfile <- "uniCox.txt"            

#read file
gene.group <- read.table(Genefile,header=T,sep="\t")
gene.exp <- read.table(GeneExpfile,header=T,sep="\t",row.names=1)
gene.exp=t(gene.exp)
gene.cox <- read.table(Coxfile,header=T,sep="\t")

#Intersect gene
colnames(gene.group) <- c('id','group')
genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))

gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

#Perpare file for network
gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)   #gene1 \t gene2 \t cor
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
gene.edge <- gene.melt[gene.melt$pvalue<0.0001,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*6

#Prepare result properties file
gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,'purple',"green")
gene.node$pvalue <- gene.cox$pvalue
# pvalue size
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- "network.node.txt"
edgefile <- "network.edge.txt"
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)


#Draw network 
node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

#output
pdf(file="network.pdf", width=10, height=8)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

#Node coordinates
coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
  if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree = degree.degree[as.numeric(degree.cut)]

#Define the pie chart, the color of the left semicircle represents the type of m6A, 
#and the right semicircle represents the gene high-risk gene or low-risk gene
values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

#draw
plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
     vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
     vertex.label.cex=2,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.4,
     vertex.color=V(g)$color,vertex.label.dist=1,vertex.label.degree=degree)

par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)
#add risk legend
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c('purple','green'),pch=16,bty="n",cex=2.5)
#add p-value
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.0001','Negative correlation with P<0.0001'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()

#m6A_cluster
library(NMF)
expFile="m6aGeneExp.txt"            
workDir="/Users/qianruzhu/Desktop/HZMUST/Fig2/m6Acluster"    
setwd(workDir)      
#read expression file
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
#NMF
res=nmf(data, rank=2:10, method="brunet", nrun=10, seed=123456)
pdf(file="cophenetic.pdf", width=8, height=7, onefile=F)
plot(res)
dev.off()
#draw
pdf(file="heatmap.all.pdf", width=15, height=15, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix",
             info=FALSE)
dev.off()
#the number of cluster
clusterNum=3       
res=nmf(data, rank=clusterNum, method="brunet", nrun=10, seed=123456)
Cluster=predict(res)
Cluster=as.data.frame(Cluster)
Cluster$Cluster=paste0("C", Cluster$Cluster)
clusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(clusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)
#heatmap
pdf(file="heatmap.pdf", width=6, height=6, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix", 
             info=FALSE)
dev.off()

#Fig2_2: m6A_clusterSurvival
library(survival)
library(survminer)
clusterFile="m6aCluster.txt"     
cliFile="clinical_804_RFS.txt"               
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig2/m6Acluster_survival")     

#read expression and clinical file
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("RFSevent", "RFS")
cli$RFS=cli$RFS/12

#data merge
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#Survival Difference Statistics
length=length(levels(factor(rt$m6Acluster)))
diff=survdiff(Surv(RFS, RFSevent) ~ m6Acluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(RFS, RFSevent) ~ m6Acluster, data = rt)
#print(surv_median(fit))

#draw
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="m6Acluster",
                   legend.labs=levels(factor(rt[,"m6Acluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="Survival_RFS.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()

#---------------------------------------------------------Fig3------------------------------------------------
#Fig3_1: m6A_cluster_GSVA

#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
expFile="merge_804.txt"            
clusterFile="cluster.txt"    
gmtFile="h.all.v2022.1.Hs.symbols.gmt"               
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig3/GSVA")

#read expression file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=t(rt)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
#GSVA analysis
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
geneSets <- getGmt("h.all.v2022.1.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#read cluster file
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#data merge
gsvaResult=read.table("gsvaOut.txt", header=T, sep="\t", check.names=F)
row.names(gsvaResult)=gsvaResult$id
gsvaResult=gsvaResult[,-1]
gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
gsvaCluster=cbind(gsvaCluster, Project)

#Difference analysis
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$m6Acluster)
comp=combn(levels(factor(allType)), 2)
for(i in 1:ncol(comp)){
  #Sample grouping
  treat=gsvaCluster[gsvaCluster$m6Acluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$m6Acluster==comp[1,i],]
  data=rbind(con, treat)
  #Difference analysis
  Type=as.vector(data$m6Acluster)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #Output the difference of all pathway
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #output significant difference
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  #diffSig=read.table("C_A_immune_diff.txt", header=T, sep="\t", check.names=F, row.names=1)
  
  #cluster color
  bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["m6Acluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  
  #Plot differential pathways heatmap
  termNum=20
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap0.pdf"),height=6,width=10)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
           cluster_cols =F,
           show_colnames = F,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 10,
           fontsize_row=7,
           fontsize_col=10)
  dev.off()
}
#Bar plot
library(tidyverse)  
library(ggthemes)
#install.packages("ggprism")
library(ggprism)
P.value=0.05

degs <- read.table("C3-C2.all.txt", header=T, sep="\t", row.names=1, check.names=F) #载入gsva的差异分析结果
Diff <- rbind(subset(degs,logFC>0)[1:20,], subset(degs,logFC<0)[1:20,]) #选择上下调前20通路     
dat_plot <- data.frame(id  = row.names(Diff),
                       p   = Diff$P.Value,
                       lgfc= Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$p <- as.numeric(dat_plot$p)
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$id <- str_replace(dat_plot$id, "KEGG_","");dat_plot$id[1:10]
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= P.value,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels=c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## Set label
low1 <- dat_plot %>% filter(lg_p < log10(P.value)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log2(P.value)) %>% nrow()
high1 <- nrow(dat_plot)

p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
                                fill = threshold)) +
  geom_col()+
  coord_flip() + 
  scale_fill_manual(values = c('Up'='#FF9933','Not'='#cccccc','Down'='#99CCFF')) +  
  geom_hline(yintercept = c(-log2(P.value),log2(P.value)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('-log10(P.Value) of GSVA score') + 
  guides(fill="none")+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black') + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 黑色标签

ggsave("GSVA_barplot_C3 vs C2.pdf",p,width = 15,height  = 15)


##Fig3_2:ssGSEA
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
expFile="merge_804.txt"           
gmtFile="immune.gmt"        
clusterFile="m6aCluster.txt"      
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig3/ssGSEA")

#read expression file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt <- t(rt)
rt=as.matrix(rt)
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(rt)

#read gmt file
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA analysis
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
# score ssGSEA and correction
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#out the score of ssGSEA
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#read cluster file
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F,row.names = 1)

#data merge
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), cluster$ID)
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind.data.frame(ssgseaScore, cluster)


#Convert data to ggplot2 input file
data=melt(scoreCluster, id.vars=c("m6Acluster"))
colnames(data)=c("m6Acluster", "Immune", "Fraction")

#boxplot
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"m6Acluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="m6Acluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="m6Acluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot.pdf", width=8, height=6.5)                          #输出图片文件
p+stat_compare_means(aes(group=m6Acluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()



#-------------------------------------------------------Fig4--------------------------------------------------------------
#Fig4_1: m6A_cluster_DEG

library(limma) 
library(VennDiagram)
expFile="merge_804.txt"          
cluFile="cluster.txt"      
adj.P.Val.Filter=0.001      
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig4/cluster_DEG")

#read expressin file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt <- t(rt)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#DEGs between the different m6A clusters
#read clusterfile
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#analysis the difference
logFCfilter=0
geneList=list()
cluster <- cluster$m6Acluster
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  #print(contrast)
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #Output the difference of all genes
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #output the difference result
  diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}

#venn plot
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#save intersect genes
interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#save the expression of intersect genes
interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)

#Fig4_2: interGene_Unicox

library(limma)
library(survival)
expFile="interGeneExp.txt"    
cliFile="clinical_804.txt"             
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig4/Unicox")

#read expression file
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data1=t(data)
rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))

#read clinical file
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #??ȡ?ٴ??ļ?
cli$RFS=cli$RFS/12

#merge clinical ab=nd expression file
cli=cli[,-(1:2)]
cli=cli[,-(3:4)]
rt=cbind(cli, data1)

#UniCox analysis
outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox????
  cox <- coxph(Surv(RFS, RFSevent) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#output
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

#output the significant results
sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uniSigGeneExp.txt", sep="\t", quote=F, col.names=F)


#Fig4_3: geneCluster
library(NMF)
expFile="uniSigGeneExp.txt"            
workDir="/Users/qianruzhu/Desktop/HZMUST/Fig4/geneCluster"    
setwd(workDir)      
#read file
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
#NMF
res=nmf(data, rank=2:10, method="brunet", nrun=10, seed=123456)
pdf(file="cophenetic.pdf", width=8, height=7, onefile=F)
plot(res)
dev.off()
#heatmap
pdf(file="heatmap.all.pdf", width=15, height=15, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix",
             info=FALSE)
dev.off()
#the number of cluster
clusterNum=2       
res=nmf(data, rank=clusterNum, method="brunet", nrun=10, seed=123456)
Cluster=predict(res)
Cluster=as.data.frame(Cluster)
Cluster$Cluster=paste0("C", Cluster$Cluster)
clusterOut=rbind(ID=colnames(Cluster), Cluster)
write.table(clusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)
#heatmap
pdf(file="heatmap.pdf", width=6, height=6, onefile=F)
consensusmap(res,
             annRow=NA,
             annCol=NA,
             #tracks=c("consensus:"),
             main="Consensus matrix", 
             info=FALSE)
dev.off()

#genecluster_survival
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig4/geneCluster")

library(survival)
library(survminer)
clusterFile="cluster.txt"   
cliFile="clinical_804.txt" 

#read file
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,-(1:2)]
cli <- cli[,-(3:4)]
colnames(cli)=c("fustat", "futime")
cli$futime=cli$futime/12

#data merge
rt=cbind(cli,cluster)

#RFS
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#RFS plot
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"Cluster"])))]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend.title="Cluster",
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival_RFS_genecluster.pdf", onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()

#read file
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli <- cli[,-(1:4)]
colnames(cli)=c("fustat", "futime")
cli$futime=cli$futime/12

#OS analysis
length=length(levels(factor(rt$Cluster)))
diff=survdiff(Surv(futime, fustat) ~ Cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Cluster, data = rt)
#print(surv_median(fit))

#OS plot
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(rt[,"Cluster"])))]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend.title="Cluster",
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="survival_OS_genecluster.pdf", onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()

#PCAscore
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig4/PCAscore")
expFile="uniSigGeneExp.txt"   
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

pca=prcomp(data, scale=TRUE)
value=predict(pca)
m6Ascore=value[,1]+value[,2]
m6Ascore=as.data.frame(m6Ascore)
scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
write.table(scoreOut, file="m6Ascore.txt", sep="\t", quote=F, col.names=F)

library(survival)

scoreFile="m6Ascore.txt"    
cliFile="clinical_804.txt"           

#merge dataframe
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
#RFS
cli <- cli[,-(1:2)]
cli <- cli[-(3:4)]
#OS
cli <- cli[-(1:4)]
colnames(cli)=c("fustat", "futime")
cli$futime=cli$futime/12
data <- cbind.data.frame(cli,score)

#cutoff
res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("m6Ascore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"m6Ascore"]<=cutoff, "Low", "High")
data$group=Type
outTab=rbind(id=colnames(data), data)
write.table(outTab, file="m6Ascore.group.txt", sep="\t", quote=F, col.names=T)

#RFS/OS survival
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
#print(surv_median(fit))

#plot
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=data,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="m6Ascore",
                   legend.labs=levels(factor(data[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=12,
                   xlab="Time(years)",
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="m6Arisk_RFS.pdf", onefile = FALSE, width=7, height=5.5)
print(surPlot)
dev.off()

#fig4_cibersort
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig4/cibersort")
library(e1071)
library(parallel)
library(preprocessCore)
library(limma)
#read GSE expression file
rt=read.table("GSE39582.txt", header=T, sep="\t", row.names=1,check.names=F)
exp=as.matrix(rt)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#data convertion
v=voom(rt, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(data,file="uniq.symbol.txt",sep="\t",quote=F,col.names=T)        #输出文件

# run CIBERSORT
source("irgCMAP24.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

#plot
library(limma)
library(pheatmap)
library(ggpubr)
library(vioplot)
#read cibersort result
immFile="CIBERSORT-Results.txt"   #免疫细胞结果浸润文件 
riskFile="m6Ascore.group.txt"          #风险文件  
pFilter=0.05                      #免疫细胞浸润结果过滤条件

#read input file
immune=read.table(immFile, header=T, sep="\t", check.names=F,row.names = 1)
data=as.matrix(immune[,1:(ncol(immune)-3)])

#read clinical file
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=rownames(risk)[risk[,"group"]=="Low"]
highSample=rownames(risk)[risk[,"group"]=="High"]

#high risk vs low risk
lowSameSample=intersect(row.names(data), lowSample)
highSameSample=intersect(row.names(data), highSample)
data=t(data[c(lowSameSample,highSameSample),])
data <- as.matrix(data)
conNum=length(lowSameSample)
treatNum=length(highSameSample)
write.table(data,file="immune.txt",sep="\t",quote=F,col.names=T)        #输出文件
highSameSample=as.data.frame(highSameSample)
lowSameSample=as.data.frame(lowSameSample)
write.table(highSameSample,file="highSameSample.txt",sep="\t",quote=F,col.names=T)  
write.table(lowSameSample,file="lowSameSample.txt",sep="\t",quote=F,col.names=T)  

#barplot
pdf("barplot.pdf",height=10,width=18)
library(RColorBrewer)
col <- colorRampPalette(brewer.pal(8,"Set1"))
#col <- topo.colors(nrow(data),alpha = 1)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col(22),yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col(22),pch=15,bty="n",cex=1.3)
dev.off()

#22immune cell boxplot
library(ggplot2)
library(reshape2)

###Constructing group comparisons and immune cell expression data
sample=read.table("m6Ascore.group.txt", header=T, sep="\t", check.names=F,row.names = 1)
immune=read.table("immune.txt", header=T, sep="\t", check.names=F, row.names = 1)
immune <-as.matrix(t(immune))
samesample <- intersect(row.names(sample),row.names(immune))
sample <- sample[samesample,]
immune <- immune[samesample,]
Immune=cbind.data.frame(sample,immune)

###Convert the data to ggplot2 input file
data=melt(Immune,id.vars=c("sample"))
colnames(data)=c( "risk", "immune","Expression")

#boxplot
p=ggboxplot(data, x="immune", y="Expression", color = "risk", 
            ylab="Fraction",
            xlab="",
            legend.title="risk",
            palette = c("red","blue"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=risk),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出pdf
pdf(file="boxplot_immune.pdf", width=7, height=5)
print(p1)
dev.off()


#———————————————————————————————————————————————————————————————————————————Fig5—————————————————————————————————————————————————————————————————————————
#Fig5_1:mRNA_extract
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)

#gencode.v22.annotation.gtf
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig5/mRNA_extract")     
gtf = rtracklayer::import("gencode.v39.annotation.gtf")
class(gtf)
gtf = as.data.frame(gtf);dim(gtf)
colnames(gtf)
table(gtf$type)
#Extract gene names from file
gtf_gene = gtf[gtf$type=="gene",]
table(gtf_gene$gene_type)
gtf_gene <- gtf[gtf$type=="gene",] %>%
  select(gene_id,gene_name,gene_type)
gtf_mRNA <- gtf_gene[gtf_gene$gene_type=="protein_coding",]

#read expression file
expfile=read.table("merge_804.txt", sep="\t",header = T,check.names = F,row.names = 1)
exp <- t(expfile)

#intersect the symbol of gtf annotation file  and the rownames of expression matrix 
same_mRNA = intersect(rownames(exp),gtf_mRNA$gene_name);length(same_mRNA)

#mRNA expression matrix

mRNA_exp <-  as.data.frame(exp[same_mRNA,])
write.table(mRNA_exp, file = "mRNA_exp.txt",sep = "\t", row.names = T)

#preWGCNA
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig5/preWGCNA")     
#acquire m6A-mRNAs 
library(limma)
corFilter=0.3            #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准

#read matrix
exp=read.table("mRNA_exp.txt", header=T, sep="\t", check.names=F)
exp=as.matrix(exp)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#read m6aGeneExp
rt1=read.table("m6aGeneExp.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
m6A=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
m6A=avereps(m6A)
m6A=m6A[rowMeans(m6A)>0.1,]

#correlation test
outTab=data.frame()
for(i in row.names(data)){
  if(sd(data[i,])>0.1){
    test=wilcox.test(data[i,])
    if(test$p.value<0.05){
      for(j in row.names(m6A)){
        x=as.numeric(data[i,])
        y=as.numeric(m6A[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(m6A=j,data=i,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(m6A=j,data=i,cor,pvalue,Regulation="negative"))
        }
      }
    }
  }
}

#output network 
write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)
#output node
mRNANode=data.frame(Node=unique(as.vector(outTab[,"data"])), Type="mRNA")
m6ANode=data.frame(Node=unique(as.vector(outTab[,"m6A"])), Type="m6A")
nodeOut=rbind(mRNANode, m6ANode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)

#out m6A-mRNA EXP
m6a_mRNA=unique(as.vector(outTab[,"data"]))
m6a_mRNAexp=data[m6a_mRNA,]
m6a_mRNAexp=rbind(ID=colnames(m6a_mRNAexp), m6a_mRNAexp)
write.table(m6a_mRNAexp,file="m6a_mRNAExp.txt",sep="\t",quote=F,col.names=F)
#|cor|>0.3
EXP_all=read.table("mRNA_exp.txt", header=T, sep="\t", check.names=F,row.names = 1)
RNA=read.table("cor_0.3.txt", header=T, sep="\t", check.names=F)
EXP_m6A_related=EXP_all[RNA$mRNA,]
EXP_m6A_related=t(EXP_m6A_related)
write.table(EXP_m6A_related,file="EXP_m6A_related.txt",sep="\t",quote=F,col.names=T)

#WGCNA  
library(limma)
library(WGCNA)
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig5/WGCNA")     
expFile="EXP_m6A_related.txt"   
cliFile="clinical_804.txt"      
   
#read expfile
exp <- read.table(expFile, header=T, sep="\t", row.names=1, check.names=F)
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=log2(data+1)
data=data[apply(data,1,sd)>0,]    
datExpr0=avereps(data)

###delete NAs
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

####  sample clustering
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###cutHeight
abline(h = 10000, col = "red")
dev.off()

###cutHeight
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]

###power
allowWGCNAThreads()    #多线程工作
powers = c(1:20)        #幂指数范围
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="2_scale_independence.pdf",width=9,height=6)
par(mfrow = c(1,2))
cex1 = 0.9
###Fitting index and power value scatterplot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 
### Mean Connectivity and Power Scatter Plot
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
###sft
sft #the best powerֵ value
softPower =sft$powerEstimate 
adjacency = adjacency(datExpr0, power = softPower)
softPower

###TOM matrix
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

###gene clustering
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="3_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


###Dynamic Clipping Module Recognition
minModuleSize=30     
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="4_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###View Similar Module Clusters
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="5_Clustering_module.pdf",width=7,height=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.3  #???и߶ȿ??޸?
abline(h=MEDissThres, col = "red")
dev.off()


###Merge similar modules
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="6_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


###Module and Trait Data Heatmap
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cli), row.names(MEs))
MEs=MEs[sameSample,]
datTraits=cli[sameSample,]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="7_Module_trait.pdf", width=6.5, height=5.2)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

### Genes in modules
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "module_all.txt",sep="\t",row.names=F,quote=F)


### Output the genes in each module
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

##FC
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig5/FC")     
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       
qvalueFilter=0.05       

#color
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("module_5.txt", header=T, sep="\t", check.names=F)     

#genename convert to entrezIDs
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO analysis
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#save result
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#Number
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#barplot
pdf(file="GO_barplot.pdf", width=9, height=10)
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#bubble plot
pdf(file="GO_bubble.pdf", width=9, height=12)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()



#-----------------------------——————————————————————-------Fig6------------------——————————————————————————————————————---------------------
#Unicox
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig6/RFS/Unicox")     
library(glmnet)
library(survival)
data1=read.table("EXP_m6A_related.txt",header = T,row.names=1, sep = "\t")
data1=as.matrix(t(data1))
data2=read.table("module_5.txt",header = T,sep = "\t")
samegene=intersect(row.names(data1),data2$Gene)
merge_data <- data1[samegene,]
merge_data=as.data.frame(t(merge_data))
cli=read.table("clinical_804.txt",header = T,row.names=1, sep = "\t")
#cli <- cli[,-(1:4)] #for OS analysis
cli <- cli[,-(1:2)]
cli <- cli[,-(3:4)]
colnames(cli)=c("fustat","futime")
Unicox_exp <- cbind.data.frame(cli,merge_data)
Unicox_exp <- na.omit(Unicox_exp)
Unicox=log2(Unicox_exp[,3:ncol(Unicox_exp)]+1)
Unicox_exp=cbind(Unicox[,1:2],Unicox_exp)
rt <- Unicox_exp
#UniCox analysis
outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox????
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#output
write.table(outTab,file="uniCox_results.txt",sep="\t",row.names=F,quote=F)

###COx
library(survminer)
library(survival)
setwd("/Users/qianruzhu/Desktop/HZMUST/Fig6/RFS/MultiCox")     

#OS/RFS
rt=read.table("EXP_m6A_related.txt", header=T, sep="\t", check.names=F, row.names=1)   
rt <- as.data.frame(t(rt))
key=read.table("key.txt", header=T, sep="\t", check.names=F)   
rt1=rt[key$Gene,]
rt1 <- as.data.frame(t(rt1))
cli=read.table("clinical_804.txt", header=T, sep="\t", check.names=F, row.names=1)   
#cli <- cli[,-(1:4)] foe OS analysis
cli <- cli[,-(1:2)]
cli <- cli[,-(3:4)]
colnames(cli) <- c("fustat","futime")
rt2 <- cbind.data.frame(cli,rt1)
mRNAmulcox <- rt2
mRNAmulcoxEXP=log2(mRNAmulcox[,3:ncol(mRNAmulcox)]+1)
mRNAmulcox=cbind(mRNAmulcox[,1:2],mRNAmulcoxEXP)

#OS/RFS COX
multiCox=coxph(Surv(futime,fustat) ~ ., data = mRNAmulcox)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#Cox forestplot
ggforest(multiCox, data = mRNAmulcox)
dev.new()

setwd("/Users/qianruzhu/Desktop/HZMUST/Fig6/gene")     

####Survival of the significant genes
library(survminer)
library(survival)
inputmRNA=read.table("EXP_m6A_related.txt", header=T, sep="\t", check.names=F, row.names=1)   
inputmRNA <-as.data.frame(inputmRNA)
cli=read.table("clinical_804.txt", header=T, sep="\t", check.names=F, row.names=1)  
cli <- cli[,-(1:2)]
colnames(cli) <- c('RFSevent','RFS','OSevent','OS')
exp <- as.data.frame(inputmRNA$ACVRL1)
row.names(exp) <- row.names(inputmRNA)
exp_key <- cbind.data.frame(cli,exp)
risk=factor(ifelse(inputmRNA$ACVRL1>median(inputmRNA$ACVRL1),"High","Low"))
#RFS
exp_key <- exp_key[,-(3:4)] #RFS
exp_key <- exp_key[,-(1:2)] #OS
kms=survfit(Surv(RFS,RFSevent)~risk,data = exp_key)
kmdffexp=survdiff(Surv(RFS,RFSevent)~risk,data=exp_key)
pValue=round(1-pchisq(kmdffexp$chisq,df=1),8)
pdf("RFS_ACVRL1.pdf")
plot(kms,lty=1.8,col=c("red","blue"),
     xlab = "survival time in year",ylab = "survival probabilities",
     main=paste("survival curve of risk score(P=",pValue,")",sep=""))
legend("bottomleft",c("High risk","low risk"),lty=1.8,col = c("red","blue"))
dev.off()
#OS
kms=survfit(Surv(OS,OSevent)~risk,data = expHSF4)
kmdffexp=survdiff(Surv(OS,OSevent)~risk,data=expHSF4)
pValue=round(1-pchisq(kmdffexp$chisq,df=1),8)
pdf("OS_ACVRL1.pdf")
plot(kms,lty=1.8,col=c("red","blue"),
     xlab = "survival time in year",ylab = "survival probabilities",
     main=paste("survival curve of risk score(P=",pValue,")",sep=""))
legend("bottomleft",c("High risk","low risk"),lty=1.8,col = c("red","blue"))
dev.off()
#---------------——————————————————————————----------------------End-----------------------————————————————————————----------------













