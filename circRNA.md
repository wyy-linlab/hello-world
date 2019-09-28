##CircRNA分析

###0. 数据获取

[Li Y *et al.,*](https://www.ncbi.nlm.nih.gov/pubmed/28794202) "CircHIPK3 sponges miR-558 to suppress heparanase expression in bladder cancer cells.", *EMBO Rep*, 2017 Sep;18(9):1646-1659

~~~
1.6G SRR5398213_1.fastq.gz
1.8G SRR5398213_2.fastq.gz
1.6G SRR5398214_1.fastq.gz
2.1G SRR5398214_2.fastq.gz
1.5G SRR5398215_1.fastq.gz
1.7G SRR5398215_2.fastq.gz
1.8G SRR5398216_1.fastq.gz
2.1G SRR5398216_2.fastq.gz
1.4G SRR5398217_1.fastq.gz
1.7G SRR5398217_2.fastq.gz
1.7G SRR5398218_1.fastq.gz
1.8G SRR5398218_2.fastq.gz
~~~

> 参考基因组
>
> 来源于gencode网站
>
> ~~~
> ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
> ~~~

### 1. 原始数据的质控、比对和过滤

>####1.1 FastQC
>
>0.fastqc.sh
>
>~~~shell
>fastqc 0.raw/*.fastq.gz -o qc_out.pre -t 60
>~~~
>
>### 1.2 trim
>
>1.trim.sh
>
>~~~shell
>for i in `ls 0.raw/*\_1.fastq.gz`
>do
>        i=${i/*\/}
>        i=${i/\_1.fastq.gz/}
>        echo $i
>        $trim_galore -q 25 --phred33 --length 25 -e 0.1 --stringency 4 --paired -o ./1.clean 0.raw/$i\_1.fastq.gz 0.raw/$i\_2.fastq.gz --path_to_cutadapt /home/software/anaconda3/bin/cutadapt &> ./logs/1.trim/$i.trim.log &
>done
>~~~
>
>~~~
>
>SRR5398213_1  3,213,297,842 bp (91.8%)
>SRR5398213_2  3,055,485,437 bp (87.3%)
>SRR5398214_1  4,025,361,512 bp (94.8%)
>SRR5398214_2  3,746,700,768 bp (88.3%)
>SRR5398215_1  3,160,265,775 bp (94.2%)
>SRR5398215_2  2,980,231,831 bp (88.8%)
>SRR5398216_1  3,973,917,300 bp (93.7%)
>SRR5398216_2  3,735,935,757 bp (88.1%)
>SRR5398217_1  3,359,372,050 bp (95.4%)
>SRR5398217_2  3,175,470,120 bp (90.2%)
>SRR5398218_1  3,551,385,997 bp (92.7%)
>SRR5398218_2  3,388,583,388 bp (88.5%)
>~~~
>
>#### 1.3 align
>
>2.align.sh
>
>~~~bash
>for i in `ls 1.clean/*_1_val_1.fq`
>do
>        i=${i/*\//}
>        i=${i/_1_val_1.fq/}
>        echo $i
>        $bwa mem -T 19 -t 70 /home/zhaojp/ref/index/GRCh37.fa 1.clean/$i\_1_val_1.fq 1.clean/$i\_2_val_2.fq > ./2.align/$i.sam 2> logs/2.align/$i.logs
>done
>~~~

### 2. 利用CIRI2进行circRNA识别

> 3.ciri.sh
>
> ~~~bash
> for i in `ls 2.align/*.sam`
> do
>         i=${i/*\//}
>         i=${i/.sam/}
>         echo $i
>         perl $ciri -I 2.align/$i.sam -O 3.cirrna/$i.circular.txt -F /home/zhaojp/ref/index/GRCh37.fa -A /home/zhaojp/ref/index/GRCh37.gtf -T 40 -G logs/3.ciri/$i.log
> done
> ~~~
>
> 利用以下两个脚本获取cirRNA的表达矩阵
>
> getCircularExpMatrix.pl
>
> ~~~perl
> use strict;
> use warnings;
> #perl getCircularExpMatrix.pl S12.circular.txt S13.circular.txt S14.circular.txt S15.circular.txt
> 
> my %hash=();
> my @fileIndexs=();
> my %geneHash=();
> 
> foreach my $file(@ARGV)
> {
> 	my $fileIndex=$file;
> 	$fileIndex=~s/(.+?)\..+/$1/g;
> 	push(@fileIndexs,$fileIndex);
> 	open(RF,"$file") or die $!;
> 	while(my $line=<RF>)
> 	{
> 		next if($.==1);
> 		next if($line=~/^\n/);
> 		my @arr=split(/\t/,$line);
> 		my @zeroArr=split(/\:|\|/,$arr[0]);
> 		$zeroArr[1]=$zeroArr[1]-1;
> 		my $circularPos="$zeroArr[0]:$zeroArr[1]-$zeroArr[2]" . $arr[10];
> 		${$hash{$circularPos}}{$fileIndex}=$arr[4];
> 		$geneHash{$circularPos}=$arr[9];
> 	}
> 	close(RF);
> }
> 
> open(WF,">circularExpMatrix.txt") or die $!;
> print WF "circRNA_ID\t" . join("\t",@fileIndexs) . "\tGene\n";
> foreach my $key(keys %hash)
> {
> 	my %valueHash=%{$hash{$key}};
> 	print WF $key;
> 	foreach my $fileIndex(@fileIndexs)
> 	{
> 		if(exists $valueHash{$fileIndex})
> 		{
> 			print WF "\t$valueHash{$fileIndex}";
> 		}
> 		else
> 		{
> 			print WF "\t0";
> 		}
> 	}
> 	print WF "\t$geneHash{$key}\n";
> }
> close(WF);
> ~~~
>
> cirBaseAnnotation.pl
>
> ~~~perl
> use strict;
> use warnings;
> #perl cirBaseAnnotation.pl /home/lexb/software/rna/circularRNA/reference/cirbase/human_hg19_circRNAs_putative_spliced_sequence.fa circularExpMatrix.txt circularAnnMatrix.txt
> 
> my %hash=();
> open(RF,"$ARGV[0]") or die $!;
> while(my $line=<RF>)
> {
> 	if($line=~/^>(.+?)\|(.+?)\|/)
> 	{
> 		if(exists $hash{$2})
> 		{
> 			print $2 . "\n";
> 		}
> 		else
> 		{
> 			$hash{$2}=$1;
> 		}
> 	}
> }
> close(RF);
> 
> open(RF,"$ARGV[1]") or die $!;
> open(WF,">$ARGV[2]") or die $!;
> while(my $line=<RF>)
> {
> 	chomp($line);
> 	if($.==1)
> 	{
> 		print WF $line . "\tcircBase_ID\n";
> 		next;
> 	}
> 	my @arr=split(/\t/,$line);
> 	if(exists $hash{$arr[0]})
> 	{
> 		print WF $line . "\t" . $hash{$arr[0]} . "\n";
> 	}
> 	else
> 	{
> 		print WF $line . "\tno\n";
> 	}
> }
> close(WF);
> close(RF);
> ~~~
>
> [circularExpMatrix.txt](./result/circularExpMatrix.txt)
>
> [circularAnnMatrix.txt](./result/circularAnnMatrix.txt)

### 3. 进行基因差异表达分析以及GO富集

利用DESeq2进行差异表达分析

**3.1 DEG.R**

~~~R
workdir="~/zhaojp/cirRNA_seq/new"
setwd(workdir)
library(DESeq2)
library(ggplot2)
library('pheatmap')
data_raw=read.table("circularAnnMatrix.txt",header = TRUE,row.names = 1,sep = "\t")
mycounts = data_raw[,1:6]
#rownames(mycounts) = data_raw[,1]
sample_info=read.table("samples.info",header = TRUE,row.names = 1)
sample_info$condition = factor(c(rep("cancer",3),rep("normal",3)))

dds <- DESeqDataSetFromMatrix(mycounts, sample_info, design= ~ condition)
dds
keep = rowSums(counts(dds)) >= 1
dds = dds[keep,]
dds <- DESeq(dds)
#res = results(dds, contrast=c("condition", "normal", "cancer"),lfcThreshold=1,alpha=0.05)
res = results(dds, contrast=c("condition", "cancer", "normal"))
summary(res)
?results
table(res$padj <0.05)

resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

res1 <- res[order(res$padj),]
resdata <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)),
                 by="row.names",sort=FALSE)
write.csv(resdata,file = "cancer_vs_normal.csv")

fdr = 0.05
logFC = 1

tT=resOrdered
tT=na.omit(tT)

tT[fdr > tT[,"padj"] & tT[,"log2FoldChange"] >= logFC, ncol(tT)+1] = "Up"
tT[fdr > tT[,"padj"] & -logFC >= tT[,"log2FoldChange"], ncol(tT)] = "Down"
tT[tT[,"padj"] >= fdr | logFC > abs(tT[,"log2FoldChange"]),ncol(tT)] = "Normal"
colnames(tT)[ncol(tT)] = "Regulate"
write.table(tT, file = "DEG.txt",  quote = F, sep = "\t")

deg = tT[fdr > tT[,"padj"]  & abs(tT[,"log2FoldChange"]) >=1,]
write.table(deg,file="deg1.txt",quote = F, sep = "\t") 

temp1 = tT[,c("padj","log2FoldChange","Regulate")]
temp1[,"padj"] = -log10(temp1$padj)
colnames(temp1)=c("-log10padj","logFC","Regulate")
temp1$Regulate=factor(temp1$Regulate, levels=c("Up","Down","Normal"), order=T)


P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10padj"]))+
  geom_point(aes(color=temp1$Regulate))+
  scale_color_manual(values =c("Up" = "red", "Down" = "blue", "Normal" = "grey"))+
  labs(x="log2FC",y="-log10FDR")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  ylim(0,15)+
  theme(plot.title = element_text(size = 25,face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),
        legend.position = 'right',
        legend.key.size=unit(0.8,'cm'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20,face = "bold", vjust = 0.5, hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = "black"))

pdf(file="volcano.pdf", width=12, height=10)
P_volcano
dev.off()


###?????
col = colorRampPalette(c('blue',  'white','red'))(256)
scal = "row" ####("none","row","columa")
annotation_col = data.frame(sample_info$condition)


rld <- rlog(dds, blind=FALSE)
eSet = assay(rld)

n.sample=ncol(eSet)
cols <- rainbow(n.sample*1.2)
boxplot(eSet, col = cols,main="expression value",las=2)
hist(eSet)

DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
deg_eset = eSet[match(DEG_list,rownames(eSet)),]


rownames(annotation_col) <- colnames(deg_eset)
ann_colors = list(Group=c(differentiated="green",undifferentiated="yellow"))

pdf(file=paste("1.deg_expressionHeatmap.pdf",sep=""), 
    width=9, height=12, onefile = FALSE)
pheatmap(deg_eset, color = col, 
         cluster_rows = T, cluster_cols=T, 
         scale = scal, show_rownames = F, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
dev.off()

vsd <- vst(dds,blind = FALSE)
pcaData = plotPCA(vsd, intgroup = c( "condition", "type"), returnData = TRUE)
#plotPCA(vsd,intgroup = c( "condition", "type"), "batch")
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("pca.pdf",width=6, height=6)
ggplot(pcaData, aes(x = PC1, y = PC2, color = type, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
dev.off()
~~~

<table><tr>
<td><img src=./result/0.samples_info.pdf border=0></td>
<td><img src=./result/pca.pdf border=0></td>
</tr></table>

<table><tr>
<td><img src=./result/1.deg_expressionHeatmap.pdf border=0></td>
<td><img src=./result/volcano.pdf border=0></td>
</tr></table>

共获取了17个上调基因，55个下调基因：

[deg.txt](./result/deg1.txt)

**GO富集分析**

~~~R
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(gridExtra)
list_gid_up = deg[which(deg$Regulate=="Up"),]
list_gid_down = deg[which(deg$Regulate == "Down"),]


gid_down = sub("[.].*","",list_gid_down$Gene)
gid_up = sub("[.].*","",list_gid_up$Gene)

db = org.Hs.eg.db
org = "hsa"
deg_list_up = bitr(gid_up,fromType="ENSEMBL",toType="ENTREZID",OrgDb=db,drop = TRUE)
deg_list_down = bitr(gid_down,fromType="ENSEMBL",toType="ENTREZID",OrgDb=db,drop = TRUE)

DEG_list = list(deg_list_up$ENTREZID,deg_list_down$ENTREZID)
names(DEG_list) = c("Up","Down")
classify = "BP"

ego = compareCluster(geneClusters=DEG_list,
                     fun="enrichGO",
                     ont = classify,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = T)
pdf(file = paste(classify,"_compareCluster_dotplot.pdf", sep = ""),width=12, height=12)
dotplot(ego,showCategory=15)
dev.off()
write.table(ego,file = paste(classify,"_GO_compareCluster.txt"),quote = F, sep = "\t")

ego_down = enrichGO(gene = gid_down,
                    keyType = 'ENSEMBL',
                    ont = "BP",
                    OrgDb = db,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = T)

write.table(ego_down,file = paste(classify,"_GO_enrichment_down.txt"),quote = F, sep = "\t")

pdf(file = paste(classify,"_ego_down.pdf",sep = ""),width=24, height=12)
p1=dotplot(ego_down,title = paste0("Enrichment Go_dot of BP"),font.size = 18)
p2=barplot(ego_down, showCategory=20, x = "GeneRatio")
grid.arrange(p1,p2,ncol=2,top="Down")
dev.off()

ego_up = enrichGO(gene = gid_up,
                  keyType = 'ENSEMBL',
                  ont = "BP",
                  OrgDb = db,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = T)

write.table(ego_up,file = paste(classify,"_GO_enrichment_up.txt"),quote = F, sep = "\t")
pdf(file = paste(classify,"_ego_up.pdf",sep = ""),width=24, height=12)
p1=dotplot(ego_up,title = paste0("Enrichment Go_dot of BP"),font.size = 18)
p2=barplot(ego_up, showCategory=20, x = "GeneRatio")
grid.arrange(p1,p2,ncol=2,top="Up")
dev.off()
~~~

><div align="center">
><img src="./result/BP_ego_up.pdf">
><img src="./result/BP_ego_Down.pdf">
><img src="./result/BP_compareCluster_dotplot.pdf">
></div>

[GO_BP_down.txt](./result/BP _GO_enrichment_down.txt)

[GO_BP_up.txt](./result/BP _GO_enrichment_up.txt)

[GO_BP_compareCluster.txt](./result/BP _GO_compareCluster.txt)

