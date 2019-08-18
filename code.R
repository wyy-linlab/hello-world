library(DESeq2)
library(ggplot2)
library(EDASeq)
library(stringr)
library(doParallel)   #  multithreading
library(variancePartition)
library(ggrepel)      #  label
library(devtools)
library(easyGgplot2)  #  Put multiple graphs on the same page    #  install_github("kassambara/easyGgplot2")
library(acPCA)
library(limma)
library(edgeR)
library(dplyr)

################################
# my function
###
plotDensity<-function(data){
  plotDa<-data.frame(ratio=c(data[gene,2]/apply(data[gene,c(1,3,4)],1,mean),data[gene,5]/apply(data[gene,c(1,3,4)],1,mean),data[gene,6]/apply(data[gene,c(1,3,4)],1,mean),data[gene,7]/apply(data[gene,c(1,3,4)],1,mean)),sample=c(rep("DS2U",length(gene)),rep("DSP",length(gene)),rep("2DS3",length(gene)),rep("DS1",length(gene))))
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),plot.title = element_text(hjust = 0.5))
  ggplot(plotDa,aes(x=ratio,color=sample))+geom_density()+theme
}

plotPca<-function(dataReg){
  pcaOut <- as.data.frame(prcomp(dataReg,scale. =T  )$rotation)
  pcaOut$type<-sample$type
  percentage <- round(prcomp(dataReg)$sdev^2 / sum(prcomp(dataReg)$sdev^2) * 100, 2)
  percentage <- paste( colnames(pcaOut), "(", paste( as.character(percentage), "%", ")", sep=""))
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),plot.title = element_text(hjust = 0.5))
  ggplot(pcaOut,aes(x=PC1,y=PC2,color=type,label=row.names(pcaOut)))+geom_point(size=3)+xlab(percentage[1]) + ylab(percentage[2])+theme+geom_text_repel(aes(x=PC1,y=PC2,label=row.names(pcaOut)))+xlim(-1,1)+ylim(-1,1)
}

countToFpkm <- function(counts,lengths){sum=apply(counts,2,sum);sweep(counts*1E9/lengths,2,sum,"/")}
fpkmToTpm <- function(fpkm){sum=apply(fpkm,2,sum);sweep(fpkm*1E6,2,sum,"/")}


################################

setwd("~/wjs/project/others/liuLab")
data<- read.table("data/data.txt",sep="\t",row.names = 1,header=T,check.names = F)[,8:14]
sample <- data.frame(row.names = colnames(data),type = c(rep("WT",4),rep("DS",3)),individual=c(1:7))

###########################
# make cibersortx data
###########################
#tang
tangSample<-rbind(read.table("data/cibersortx/tang/PFC_astrocyte_rf_ident.txt",sep="\t",row.names = 1,header=T),rbind(read.table("data/cibersortx/tang/PFC_ExcitatoryNeurons_rf_c7.txt",sep="\t",row.names = 1,header=T), rbind(read.table("data/cibersortx/tang/PFC_interneurons_rf_c8_monocle.txt",sep="\t",row.names = 1,header=T)[,c(3,4)],  rbind(read.table("data/cibersortx/tang/PFC_microglia_rf_ident.txt",sep="\t",row.names = 1,header=T),  rbind(read.table("data/cibersortx/tang/PFC_NPCs_rf_c9_monocle_ident.txt",sep="\t",row.names = 1,header=T)[,c(3,4)], read.table("data/cibersortx/tang/PFC_OPC_rf_c4_ident.txt",sep="\t",row.names = 1,header=T))))))
tangSample[,1]<-c(rep(c("astrocyte","ExcitatoryNeurons","interneurons","microglia","NPC","OPC"),c(76,1055,701,68,288,117)))
#setdiff(row.names(sample),c(tangSample))   #GW09_PFC1_sc3  Stem cells GW09 #GW10_PFC3_sc46 Stem cells GW10 # GW16_PFC1_DL3_sc14 Neurons GW16 # GW16_PFC1_L1_sc1      Neurons GW16   
tang<-t(read.table("data/cibersortx/GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls",sep="\t"))
tmp <- t(merge(tang,tangSample,by=0))  ; colnames(tmp)<-tmp[1,] ; tmp[1,]<-tmp[24155,]  ;  tang<-tmp[-c(24155,24156),]
write.table(tang,"data/cibersortx/processData/tangTpm.txt",sep="\t",col.names = F,quote=F)
# DSCAM
tmp<-tang[c(1,which(row.names(tang)=="TCF4")),] ;row.names(tmp)[1]<-"type"
tmp1<-merge(t(tmp),tangSample,by=0)
#tmp[2,as.numeric(tmp[2,])>50]<-NA;tmp<-t(na.omit(t(tmp)))
#ylims <- data.frame(tpm=as.numeric(as.matrix(tmp)[2,]),type=as.character(as.matrix(tmp)[1,])) %>%group_by(type) %>%summarise(Q1 = quantile(as.numeric(as.matrix(tmp)[2,]), 1/4), Q3 = quantile(as.numeric(as.matrix(tmp)[2,]), 3/4)) %>%ungroup() %>%summarise(lowQ1 = min(Q1), highQ3 = max(Q3)) 
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5,colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),plot.title = element_text(hjust = 0.5))
#ggplot(data.frame(log2Tpm=log2(as.numeric(as.matrix(tmp)[2,])+1),type=as.character(as.matrix(tmp)[1,])),aes(type,log2Tpm,color=type))+geom_violin(trim = F)+geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.2)+ggtitle("DSCAM in PMID:29539641")+theme
#ggplot(data.frame(tpm=as.numeric(as.matrix(tmp)[2,]),type=as.character(as.matrix(tmp)[1,])),aes(type,tpm,color=type))+geom_violin(trim = F)+ggtitle("DSCAM in PMID:29539641")+theme+geom_boxplot(width=0.1,outlier.colour = NA)#+coord_cartesian(ylim=as.numeric(ylims)*1.05)
#ggplot(data.frame(tpm=as.numeric(as.matrix(tmp1)[,3]),timePoint=paste(as.character(as.matrix(tmp1)[,2]),as.character(as.matrix(tmp1[,5]),sep="_")),type=as.character(as.matrix(tmp1)[,2])),aes(timePoint,tpm,color=type))+ggtitle("DSCAM in PMID:29539641")+theme+geom_boxplot(outlier.shape = NA)+ylim(1,100)
ggplot(data.frame(tpm=as.numeric(as.matrix(tmp1)[,3]),timePoint=paste(as.character(as.matrix(tmp1)[,2]),as.character(as.matrix(tmp1[,5]),sep="_")),type=as.character(as.matrix(tmp1)[,2])),aes(timePoint,tpm,color=type))+ggtitle("DSCAM in PMID:29539641")+theme+geom_boxplot()


#liu count-tpm
liuCount <- read.table("data/data.txt",sep = "\t",row.names = 1,header = T,check.names = F)[,c(8:15)]
liuCount[,1:7]<-round(liuCount[,1:7])   ;  liuCount[,8]<-as.character(liuCount[,8])
dds<- DESeqDataSetFromMatrix(liuCount[,1:7],sample,~type)
dds <- estimateSizeFactors(dds)
liuCount[,1:7] <- counts(dds,normalized=TRUE)

#geneLengthHg19<-getGeneLengthAndGCContent(as.character(row.names(liuCount)),"hg19","org.db")
#geneLength<-na.omit(geneLengthHg19)
liuCount2<-liuCount[row.names(geneLength),]
liuCount2<-liuCount2[!liuCount2[,8]%in%unique(liuCount2[duplicated(liuCount2[,8]),8]),]

lengths=as.numeric(as.character(geneLength[row.names(liuCount2),1]))
liuFpkm<-countToFpkm(liuCount2[,1:7],lengths)
liuTpm <-cbind(liuCount2[,8],fpkmToTpm(liuFpkm))
colnames(liuTpm)<-c("gene","IMR90","DS2U","NC1","H9","DSP","2DS3","DS1")
write.table(liuTpm,"data/cibersortx/processData/liuTpm.txt",row.names = F,sep="\t",quote=F)

##################################
# run cibersortx
##################################  

##################################
# variancePartition
##################################
cc<-read.table("~/wjs/project/others/liuLab/result/2019-06-23/cibersortx/tang.csv",sep=",",header=T,row.names = 1)
sample1<-sample  ; sample1$NPC <- cc$NPC ; sample1$ExcitatoryNeurons <-cc$ExcitatoryNeurons;  sample1$interneurons <- cc$interneurons ; sample1$astrocyte <- cc$astrocyte
row.names(sample1)<-colnames(liuTpm)[-1]
tmp<-liuTpm[,-1]
form <- ~ type +ExcitatoryNeurons + NPC  +astrocyte
varPart <- fitExtractVarPartModel(tmp[rowSums(tmp>1) >= 0.5*ncol(tmp),], form, sample1)
plotVarPart(sortCols(varPart))

#################################
#ACPCA
#################################
tpm<-liuTpm ; row.names(tpm)<-tpm[,1]; tpm <- tpm[,-1]
acplData<-apply(tpm,1,as.numeric)
row.names(acplData)<-colnames(tpm)
resultTune <- acPCAtuneLambda(X=acplData, Y=as.matrix(sample1[,5:6]), nPC=2, lambdas=seq(0,200, 0.05),anov=F, kernel = "linear",perc = 0.05, quiet=T)
result <- acPCA(X=acplData, Y=as.matrix(sample1[,5:6]), lambda=resultTune$best_lambda,kernel="linear", nPC=2)
acPcaOut <- as.data.frame(result$Xv)
colnames(acPcaOut)<-c("PC1","PC2")
acPcaOut$type<-sample1$type
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),plot.title = element_text(hjust = 0.5))
ggplot(acPcaOut,aes(x=PC1,y=PC2,color=type))+geom_point(size=3)+theme+ggtitle("AC-PCA")+geom_text_repel(aes(x=PC1,y=PC2,label=row.names(acPcaOut)))+xlim(-6000,6000)+ylim(-5000,5000)


################################# 
#linear model
#################################
tpm<-liuTpm ; row.names(tpm)<-tpm[,1]; tpm <- tpm[,-1]
type <- 2-as.numeric(as.factor(sample1$type)) ; NPC <- as.numeric(sample1$NPC) ; interneurons <- as.numeric(sample1$interneurons) ; astrocyte <- as.numeric(sample1$astrocyte)
regvars1 <- as.data.frame(cbind(type,astrocyte))
tpmReg1 <- matrix(NA,nrow=nrow(tpm),ncol=ncol(tpm))
rownames(tpmReg1) <- rownames(tpm)   ;  colnames(tpmReg1) <- colnames(tpm)
coefmat1 <- matrix(NA,nrow=nrow(tpm),ncol=ncol(regvars1)+1)
for (i in 1:nrow(tpm)) {
  lmmod1 <- lm(as.numeric(tpm[i,])~astrocyte+type,data=regvars1)
  coef1 <- coef(lmmod1)
  coefmat1[i,] <- coef1
  tpmReg1[i,] <- as.numeric(tpm[i,] - coef1["astrocyte"]*regvars1[,"astrocyte"])
}
#barplot
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),plot.title = element_text(hjust = 0.5))
ggplot(data.frame(tpm=round(c(as.numeric(tpm["DSCAM",]),as.numeric(tpmReg1["DSCAM",])),2),sample=c(colnames(tpm),colnames(tpm)),group=c(rep(c("raw","lm"),c(7,7)))),aes(x=sample,y=tpm,fill=group))+theme+ggtitle("DSCAM_lm-raw_tpm")+geom_bar(stat = "identity",position=position_dodge())+geom_text(aes(label=tpm),vjust=1.6,color="red",position = position_dodge(0.9),size=3.5)+scale_fill_brewer(palette = "Blues")

##chr21 gene density
hg19chr21<-read.table("~/wjs/reference/human/hg19/gtf/hg19chr21.txt",sep=" ")
hg19chr21[,3]<-word(hg19chr21[,3],1,sep=fixed("."))   ;   
hg19chr21[,4]<-word(hg19chr21[,4],1,sep=fixed(";"))
row.names(hg19chr21)<-paste(hg19chr21[,3],hg19chr21[,4],sep="|")
gene<-intersect(hg19chr21[,4],row.names(tpmReg1))
plotDensity(tpmReg1)+ggtitle("adjust astrocyte chr21 log2(DS/mean(DT)) density")+xlim(-1,5)

#################################
##limma 
#################################
tpm<-liuTpm ; row.names(tpm)<-tpm[,1]; tpm <- tpm[,-1]
design <- model.matrix(~0+factor(sample1$type,levels = c("DS","WT")) +sample1$astrocyte)
rownames(design) <- colnames(tpm) ; colnames(design) <- c("DS","WT","astrocyte")
fit <- lmFit(tpm, design)
cont.matrix=makeContrasts("DS-WT",levels = design)
fit <- contrasts.fit(fit,cont.matrix)
fit <- eBayes(fit)
output1 <- topTable(fit,coef=1,n=Inf)
DEG<-na.omit(output1)
head(DEG)



