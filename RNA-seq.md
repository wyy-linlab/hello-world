##RNA-seq分析
###准备

1. 参考基因组

	人类基因组(gencode数据库)--transcripts.fa 用于建库
	   	
	~~~shell
		wget -bc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human  \
   			/release_31/gencode.v31.transcripts.fa.gz
	~~~

   获取tx2gene.info信息

2. 软件：kalliso、DESeq2
3. SRA数据或者fastq数据、样本的基本信息

	`samples.txt`

###上游分析
1. 首先利用fastq-dump --split-3拆分SRA数据
2. 利用fastqc查看fastq的数据质量 

  00.fastq.sh

 ~~~shell
  fastqc 00.data/*.fq -t 4 -o qc_out.pre
  multiqc ./qc_out.pre/
 ~~~

3. 利用trimmomatic或者trim_galore进行数据过滤获得clean数据  

  01.trim.sh
  > trimmomatic

  ~~~shell
  	for i in `ls 00.data/*_1.fastq.gz`
  	do
  		i=${i/*\//}
  		i=${i/_1.fastq.gz/}
  		echo $i
  		java -jar trimmomatic-0.39.jar PE -threads 18 -phred33 \
  			00.data/$i\_1.fastq.gz 00.data/$i\_2.fastq.gz \
  			01.clean/$i\_1.clean.fq.gz 01.clean/$i\_1.unpaired.fq.gz \
  			01.clean/$i\_2.clean.fq.gz 01.clean/$i\_2.unpaired.fq.gz \
  			ILLUMINACLIP:`PATH`/adapters/TruSeq3-PE-2.fa:2:30:10 \
  			LEADING:3 TRAILING:3 \
  			SLIDINGWINDOW:4:15 MINLEN:36 2>>./logs/log.trim
  	done	
  ~~~
  
   > trim_galore
  
    ```shell
    	for i in `ls 00.data/*_1.fastq.gz`
    	do
    		i=${i/*\//}
    		i=${i/_1.fastq.gz/}
    		echo $i
    		trim_galore -q 20 --phred33 --stringency 3 --length 20 \
    			--paired 00.data/$i\_1.fastq.gz 00.data/$i\_2.fastq.gz \
    			-o 01.data &> ./logs/log.trim.galore
    	done
    ```

4. kallisto定量处理   

  02.kallisto.sh
  >index

  ```shell
  	kallisto index -i hg38.idx *.fa
  ```
  
  > quant 
  >
  > ~~~bash
  > for i in `01.clean/*_1.clean.fastq.gz`
  >   	do
  >   		i=${i/*\//}
  >   		i=${i/_1.clean.fastq.gz/}
  >   		echo $i
  >   		kallisto quant -i ./ref/hg38.idx \
  >   		-o ./result/$i \
  >   		-t 10 -b 100 \
  >   		01.clean/$i\_1.clean.fastq.gz \
  >   		01.clean/$i\_2.clean.fastq.gz &> logs/quant.log.$i
  >   	done			
  > ~~~

5. DESeq2差异分析

   ~~~R
   	# if (!requireNamespace("BiocManager", quietly = TRUE))
   	#   install.packages("BiocManager")
   	# BiocManager::install("DESeq2")
   	# BiocManager::install("tximport")
   	# BiocManager::install("org.Hs.eg.db")
   	# BiocManager::install("reshape2")
   	# BiocManager::install("readr")
   	# BiocManager::install("rhdf5")
   	# BiocManager::install("pheatmap")
   	# BiocManager::install("reshape2")
   	# BiocManager::install("stringr")
   	# BiocManager::install("gridExtra")
   
   
     rm(list=ls())
     ### 提前安装对应的R包并正常加载
   	ibrary("DESeq2")
     library("tximport")
     library("readr")
     library('ggplot2')
     library('reshape2')
     library('pheatmap')
     
     ###设置工作路径	
     workdir="~/RNA_seq/result"
     setwd(workdir)
     	
     ###读取tx2gene、样本分组等信息
     tx2genes = read.table("tx2gene.info",sep="\t",header = TRUE)
   		#====================tx2gene.info===========================
     	#TXNAME	GENEID
     	#ENST00000456328.2	ENSG00000223972.5
     	#ENST00000450305.2	ENSG00000223972.5
     	#ENST00000488147.1	ENSG00000227232.5
     	#ENST00000619216.1	ENSG00000278267.1
     	#===========================================================
     	
     samples = read.table("samples.txt",header=T)
     	
     	#=====================samples.txt===========================
     	#  Run	       SampleName	 cell	   dex
     	#SRR1039508	GSM1275862	N61311	  untrt
     	#SRR1039509	GSM1275863	N61311	  trt
     	#SRR1039512	GSM1275866	N052611  untrt
     	#SRR1039513	GSM1275867	N052611  trt
     	#SRR1039516	GSM1275870	N080611  untrt
     	#SRR1039517	GSM1275871	N080611  trt
     	#SRR1039520	GSM1275874	N061011  untrt
     	#SRR1039521	GSM1275875	N061011  trt
     	#===========================================================
     	
     ###读取kalliso数据
     files = file.path(".",samples$Run,"abundance.tsv")
     names(files) = samples$Run
     txi = tximport(files,type="kallisto",tx2gene = tx2genes,
                  ignoreAfterBar=TRUE)
   	clodata = samples[,c("cell","dex")]
     dds = DESeqDataSetFromTximport(txi,
                                  colData = clodata,
                                  design = ~ dex)
     #head(dds)
     	
     # Prefiltering the dataset	
     	keep = rowSums(counts(dds)) >= 1
     	dds = dds[keep,]
     	
     	dds = DESeq(dds)	
     	res <- results(dds, contrast=c("dex","trt","untrt"))
     #head(res)
     	summary(res)
     	
     ###按padj的值进行排序
     	resOrdered <- res[order(res$padj),]
     	resOrdered=as.data.frame(resOrdered)
     	
     ###设定阈值
     	fdr = 0.05
     	logFC = 1
     	
     ###去除数据中的NA值
     	tT=resOrdered
     	tT=na.omit(tT)
   	###增加基因的上下调信息
     	tT[fdr > tT[,"padj"] & tT[,"log2FoldChange"] >= logFC, ncol(tT)+1] = "Up"
     	tT[fdr > tT[,"padj"] & -logFC >= tT[,"log2FoldChange"], ncol(tT)] = "Down"
     	tT[tT[,"padj"] >= fdr | logFC > abs(tT[,"log2FoldChange"]),ncol(tT)] = "Normal"
     	colnames(tT)[ncol(tT)] = "Regulate"
     	write.table(tT, file = "DEG.txt",  quote = F, sep = "\t")
     	
    	###获取上下调基因
     	deg = tT[fdr > tT[,"padj"]  & abs(tT[,"log2FoldChange"]) >=1,]
     	write.table(deg,file="deg1.txt",quote = F, sep = "\t")		
     ###根据分组差异比较结果绘制火山图
     #获得画图数据padj值和log2Fc以及上下调信息
     	temp1 = tT[,c("pvalue","log2FoldChange","Regulate")]
     	temp1[,"pvalue"] = -log10(temp1$pvalue)
     	colnames(temp1)=c("-log10pvalue","logFC","Regulate")
     	temp1$Regulate=factor(temp1$Regulate, levels=c("Up","Down","Normal"), order=T)
     #绘图
     	P_volcano=ggplot(temp1,aes(x=temp1$logFC,y=temp1[,"-log10pvalue"]))+
     			geom_point(aes(color=temp1$Regulate))+
     			scale_color_manual(values =c("Up" = "red", "Down" = "blue", "Normal" = "grey"))+
      			labs(x="log2FC",y="-log10FDR")+
      			geom_hline(yintercept=-log10(fdr),linetype=4)+
      			geom_vline(xintercept=c(-1,1),linetype=4)+
      			xlim(-5,5)+
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
     	print(P_volcano)
     	dev.off()
   ~~~

     <div align="center">
     <img src="./result/volcano.png" width = 80% height = 80% >
     </div>

~~~R
		###获得差异基因id,获取基因在样本中的表达量并绘制聚类热图
  	#加载聚类热图的R包，并设置个别参数
  	col = colorRampPalette(c('blue',  'white','red'))(256)
  	scal = "row" ####("none","row","columa")
  	annotation_col = data.frame(samples$dex)

		rld <- rlog(dds, blind=FALSE)
		eSet = assay(rld)

  	rownames(annotation_col) <- colnames(deg_eset)
  	ann_colors = list(Group=c(differentiated="green",undifferentiated="yellow"))
  	
  	#获取差异基因GeneID 
  	DEG_list = rownames(deg)
  	DEG_list = DEG_list[order(DEG_list)]
  
  	#绘制差异基因在比较组中的聚类热图
  	deg_eset = eSet[match(DEG_list,rownames(eSet)),]
  	
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
~~~

<div align="center">
<img src="./result/1.deg_expressionHeatmap.png" width =70% height = 70% >
</div>

~~~R
		###GO富集分析
  	
  	library(clusterProfiler)
  	library(org.Hs.eg.db)
  	library(stringr)
  	library(gridExtra)
  	
  	#读取文件中的GeneID (ENSEMBL) 
		#deg_file = "deg1.txt"
  	#degs = read.table(deg_file,header=T,comment.char = "",check.names=F, row.names=1)
  	list_gid_up = deg[which(deg$Regulate=="Up"),]
  	list_gid_down = deg[which(deg$Regulate == "Down"),]
  	gid_down = substring(rownames(list_gid_down),1,15)
  	gid_up = substring(rownames(list_gid_up),1,15)		
  	# GO富集分析,分"BP","CC","MF"三类进行
  	db = org.Hs.eg.db
  	org = "hsa"	
  	classify = "BP"
  	
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
  
 		 write.table(ego_up,file = paste(classify,"_GO_enrichment_down.txt"),quote = F, sep = "\t")
 		 pdf(file = paste(classify,"_ego_up.pdf",sep = ""),width=24, height=12)
 		 p1=dotplot(ego_up,title = paste0("Enrichment Go_dot of BP"),font.size = 18)
 		 p2=barplot(ego_up, showCategory=20, x = "GeneRatio")
 		 grid.arrange(p1,p2,ncol=2,top="Up")
 		 dev.off()
  deg_list_up = bitr(gid_up,fromType="ENSEMBL",toType="ENTREZID",OrgDb=db,drop = TRUE)
  deg_list_down = bitr(gid_down,fromType="ENSEMBL",toType="ENTREZID",OrgDb=db,drop = TRUE)
  
  DEG_list = list(deg_list_up$ENTREZID,deg_list_down$ENTREZID)
  DEG_list = list(gid_up,gid_down)
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
  pdf(file = paste(classify,"_compareCluster_dotplot.pdf", sep = ""),width=12, 	height=12)
  dotplot(ego,showCategory=15)
  dev.off()
	#PCA分析
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
  #pcaData
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pdf("pca.pdf",width=6, height=6)
  ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  				geom_point(size =3) +
  				xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  				ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  				coord_fixed()
  dev.off()
~~~

<div align="center">
<img src="./result/BP_ego_down.png" >
<img src="./result/BP_ego_up.png" >
<img src="./result/BP_compareCluster_dotplot.png">
<img src="./result/pca.png">
</div>

