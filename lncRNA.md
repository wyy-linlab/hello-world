## lncRNA 分析流程

### 流程图及其参考文献

<div align="center">
<img src="./result/1.流程图.png">
</div>

Huang L, Dong H, Zhou D, et al. Systematic identification of long non‐coding RNA s during pollen development and fertilization in Brassica rapa[J]. The Plant Journal, 2018, 96(1): 203-222.

### 0. 数据获取

Barr J A, Hayes K E, Brownmiller T, et al. Long non-coding RNA FAM83H-AS1 is regulated by human papillomavirus 16 E6 independently of p53 in cervical cancer cells[J]. Scientific reports, 2019, 9(1): 3662.

~~~
1.4G SRR7262883_1.fastq.gz
1.5G SRR7262883_2.fastq.gz
1.3G SRR7262884_1.fastq.gz
1.4G SRR7262884_2.fastq.gz
1.6G SRR7262885_1.fastq.gz
1.7G SRR7262885_2.fastq.gz
1.3G SRR7262886_1.fastq.gz
1.4G SRR7262886_2.fastq.gz
1.3G SRR7262887_1.fastq.gz
1.4G SRR7262887_2.fastq.gz
1.7G SRR7262888_1.fastq.gz
1.8G SRR7262888_2.fastq.gz
~~~

> 参考基因组
>
> 来源于gencode网站
>
> ~~~
> ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
> ~~~

###1. 原始数据的质控、比对和过滤

> #### 1.1 FastQC
>
> 0.fastqc.sh
>
> ~~~shell
> fastqc 0.raw/*.fastq.gz -o qc_out.pre -t 60
> ~~~
>
> #### 1.2 trim
>
> 1.trim.sh
>
> ~~~shell
> for i in `ls 0.raw/*\_1.fastq.gz`
> do
>         i=${i/*\/}
>         i=${i/\_1.fastq.gz/}
>         echo $i
>         $trim_galore -q 25 --phred33 --length 25 -e 0.1 --stringency 4 --paired -o ./1.clean 0.raw/$i\_1.fastq.gz 0.raw/$i\_2.fastq.gz --path_to_cutadapt /home/software/anaconda3/bin/cutadapt &> ./logs/1.trim/$i.trim.log &
> done
> ~~~
>
> ~~~
> 1.4G SRR7262883_1_val_1.fq.gz
> 1.5G SRR7262883_2_val_2.fq.gz
> 1.3G SRR7262884_1_val_1.fq.gz
> 1.3G SRR7262884_2_val_2.fq.gz
> 1.5G SRR7262885_1_val_1.fq.gz
> 1.6G SRR7262885_2_val_2.fq.gz
> 1.2G SRR7262886_1_val_1.fq.gz
> 1.3G SRR7262886_2_val_2.fq.gz
> 1.3G SRR7262887_1_val_1.fq.gz
> 1.3G SRR7262887_2_val_2.fq.gz
> 1.6G SRR7262888_1_val_1.fq.gz
> 1.7G SRR7262888_2_val_2.fq.gz
> SRR7262883_1  2,666,505,653 bp (95.6%)
> SRR7262883_2  2,636,723,370 bp (94.5%)
> SRR7262884_1  2,465,218,446 bp (95.2%)
> SRR7262884_2  2,436,606,642 bp (94.1%)
> SRR7262885_1  2,950,958,950 bp (95.4%)
> SRR7262885_2  2,902,368,906 bp (93.9%)
> SRR7262886_1  2,378,617,235 bp (95.2%)
> SRR7262886_2  2,356,451,505 bp (94.3%)
> SRR7262887_1  2,433,104,160 bp (95.3%)
> SRR7262887_2  2,405,935,521 bp (94.2%)
> SRR7262888_1  3,232,484,368 bp (91.9%)
> SRR7262888_2  3,193,681,971 bp (90.8%)
> ~~~
>
> #### 1.3 align
>
> 2.bulid.sh
>
> ~~~
> $path/extract_exons.py $path_ref/GRCh37.gtf > ./ref/GRCh37.exon
> $path/extract_splice_sites.py $path_ref/GRCh37.gtf > ./ref/GRCh37.ss
> cd ref 
> $path/hisat2-build -p 40 --ss GRCh37.ss --exon GRCh37.exon $path_ref/GRCh37.fa GRCh37
> ~~~
>
> 2.align.sh
>
> ~~~
> for i in `ls 1.clean/*_1_val_1.fq.gz`
> do
>         i=${i/*\//}
>         i=${i/_1_val_1.fq.gz/}
>         echo $i
>         $path/hisat2 -p 60 --dta --phred33 -x ./ref/GRCh37 -1 1.clean/$i\_1_val_1.fq.gz -2 1.clean/$i\_2_val_2.fq.gz -S 2.align/$i.sam &> logs/2.align/$i.log 
>         $samtools sort -@ 60 -o 2.align/$i.bam 2.align/$i.sam &> logs/2.align/$i.samtools.log
> done
> ~~~
>
> ~~~
> 比对率
> SRR7262883 87.99%
> SRR7262884 83.22%
> SRR7262885 83.05%
> SRR7262886 86.05%
> SRR7262887 82.94%
> SRR7262888 69.16%
> ~~~

### 2. 组装与筛选

> ####2.1 stringtie 组装
>
> 3.stringtie.sh
>
> ~~~shell
> for i in `ls 2.align/*.bam`
> do
>         i=${i/*\//}
>         i=${i/.bam/}
>         echo $i
>         $path/stringtie 2.align/$i.bam -G $path_ref/GRCh37.gtf -o 3.stringtie/$i.gtf -p 60 &> logs/3.stringtie/$i.log
> done
> ~~~
>
> ####2.2 gffcompare 筛选
>
> 4.gtffcompare.sh
>
> ~~~shell
> cd 3.stringtie
> $path/stringtie --merge -G $path_ref/GRCh37.gtf -o stringtie_merge.gtf mergelist.txt
> $path_gff/gffcompare -r $path_ref/GRCh37.gtf -o gffcmp stringtie_merge.gtf
> #awk '{print $4}' gffcmp.tracking |sort |uniq -c >statistics
> awk '{if($10 >= 200 && $6 >= 2) print $0}' gffcmp.stringtie_merge.gtf.tmap > filter_lnc.tmp
> awk '{if ($3=="u" || $3=="x" || $3=="i" || $3=="j" || $3=="o"){print $0}}' filter_lnc.tmp > filter_lnc.class
> ~~~
>
> #### 2.3 TransDecoder筛选
>
> ~~~shell
> $trans_path/util/gtf_to_alignment_gff3.pl stringtie_merge1.gtf > stringtie_merged.gff3
> $trans_path/util/gtf_genome_to_cdna_fasta.pl stringtie_merge1.gtf $ref_path/GRCh37.fa > transcripts.fasta
> 
> echo "## Extract the long ORFs"
> $trans_path/TransDecoder.LongOrfs -t transcripts.fasta
> 
> echo "## run blast"
> 
> $blast_path/makeblastdb -in ./uniprot/uniprot_sprot.fasta -dbtype prot -parse_seqids
> $blast_path/blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db ./uniprot/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 50 > blastp.outfmt6
> $blast_path/blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db ./uniprot/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 50 -soft_masking true > blastp.outfmt6.new
> 
> echo "## run pfam"
> $hmm_path/hmmpress -f ./pfam/Pfam-A.hmm
> $hmm_path/hmmscan --cpu 30 --domtblout pfam.domtblout ./pfam/Pfam-A.hmm transcripts.fasta.transdecoder_dir/longest_orfs.pep
> 
> echo "use pfam and blast results:"
> $trans_path/TransDecoder.Predict  -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
> 
> echo "## convert to genome coordinates"
> $trans_path/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 stringtie_merged.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
> ~~~

#### 3. 识别novel lncRNA

**3.1 根据transcripts.fasta.transdecoder.genome.gff3以及filter_lnc.class结果编写脚本从stringtie_merge.gtf中提取novel.gtf**

5.filter.pl

~~~perl
#!/usr/bin/perl
use strict;
use warnings;
#perl 5.filter.pl filter_lnc.class transcripts.fasta.transdecoder.genome.gff3 stringtie_merge.gtf
my ($file_class,$file_gff3,$file_raw) = @ARGV;
my (%classg,%classt,%gff3);
&get_class($file_class,\%classg,\%classt);
&get_gff3($file_gff3,\%gff3);

open OUT,">novel.gtf" or die;
open IN,"$file_raw" or die;
while (my $line = <IN>){
        chomp $line;
        next if ($line =~ /^#/);
        my $note = (split /\t/,$line)[8];
        $note =~ /gene_id\s+\"(.*?)\";\s+transcript_id\s+\"(.*?)\";/;
        my $gid = $1; my $tid = $2;
        next if exists $gff3{$gid};
        if (exists $classg{$gid} and exists $classt{$tid}){
                print OUT "$line\n";
        }
}
close IN;
close OUT;
sub get_gff3{
        my ($file,$hash) = @_;
        open IN,"$file" or die;
        while (my $line = <IN>){
                chomp $line;
                next if (!$line);
                my ($type,$note) =  (split /\t/,$line)[2,8];
                next if ($type ne "gene");
                $note =~ /ID=(.*?);/;
                my $gid = $1;
                $$hash{$gid} ++;
        #       print "$gid\n";
        }
        close IN;
}
sub get_class{
        my ($file,$hashg,$hasht) = @_;
        open IN,"$file" or die;
        while (my $line = <IN>){
                chomp $line;
                my ($gid,$tid) = (split /\t/,$line)[3,4];
                $$hashg{$gid} ++;
                $$hasht{$tid} ++;
        }
        close IN;
}
~~~

##### 3.2 CNCI 获取编码能力

6.cnci.sh

~~~shell
faToTwoBit $ref/GRCH37.fa GRCH37.2bit
python /home/software/CNCI/CNCI.py -f ./novel.gtf -g -m ve -o ./result -d ./GRCH37.2bit -p 40 &> log.cnci
python /home/software/CNCI/filter_novel_lincRNA.py -i ./result/CNCI.index -g novel.gtf -o ./result
~~~

>This step classifies the potentially novel genes into four categories according to the CNCI output: (1) novel coding genes; (2) novel lincRNA genes; (3) ambiguous genes which give rise to both coding and noncoding transcripts; (4) filtered out noncoding genes

**3.3 CPC2获取编码能力**

7.cpc.sh

~~~shell
cd 5.novel/cpc
##根据novel.gtf获取转录本序列
$trans_path/util/gtf_genome_to_cdna_fasta.pl stringtie_merge1.gtf $ref_path/GRCh37.fa > transcripts.fasta 
##运行cpc2软件
python $path/CPC2.py -i transcripts.fasta -o cpc2output.txt
~~~

**3.4 根据3.2以及3.3的结果获取novel lncRNA的并集

6.get_lncRNA_end.pl

~~~perl
#!/usr/bin/perl
use strict;
use warnings;
#print "perl 6.get_lncRNA_end.pl cpc.file cnci.file novel.gtf stringtie_merge.gtf";

my ($file_cpc,$file_cnci,$file_gtf,$file) = @ARGV;
my (%cpc,%cnci,%gtf);
&get_cpc($file_cpc,\%cpc);
&get_cnci($file_cnci,\%cnci);
&get_gtf($file_gtf,\%gtf);

my (%gid2tid,%tid2gid);

foreach my $gid (sort keys %gtf){
	my $tid = $gtf{$gid};
	if (exists $cpc{$tid} and exists $cnci{$gid}){
		$gid2tid{$gid} = $tid;
		$tid2gid{$tid} = $gid
	}
}

my $num = keys %gid2tid;
print "$num\n";

open OUT,">novel_lncRNA_end.gtf" or die;
open IN,"$file" or die;
while (my $line = <IN>){
	chomp $line;
	next if ($line =~ /^#/);
	my $note = (split /\t/,$line)[8];
	$note =~ /gene_id\s+\"(.*?)\";\s+transcript_id\s+\"(.*?)\";/;
	my $gid = $1; my $tid = $2;
	if (exists $gid2tid{$gid} and exists $tid2gid{$tid}){
		print OUT "$line\n";
	}
}
close IN;
close OUT;
sub get_gtf{
	my ($file,$hash) = @_;
	open IN,"$file" or die;
	while (my $line = <IN>){
		chomp $line;
		my ($type,$note) = (split /\t/,$line)[2,8];
		next if ($type ne "transcript");
		$note =~ /gene_id\s+\"(.*?)\";\s+transcript_id\s+\"(.*?)\";/;
		my $gid = $1; my $tid = $2;
		$$hash{$gid} = $tid;
	}
	close IN;	
}
sub get_cnci{
	my ($file,$hash) = @_;
	open IN,"$file" or die;
	while (my $line = <IN>){
		chomp $line;
		next if ($line =~ /gene_id/);
		my ($gid,$type) = (split /\t/,$line);
		if ($type eq "novel\_lincRNA"){
			$$hash{$gid} ++;
		}
	}
	close IN;
}
sub get_cpc{
	my ($file,$hash) = @_;
	open IN,"$file" or die;
	while (my $line = <IN>){
		chomp $line;
		next if ($line =~ /^#/);
		my ($tid,$type) = (split /\t/,$line)[0,-1];
		if ($type eq "noncoding"){
			$$hash{$tid} ++;
		}
	}
	close IN;
}
~~~

最后筛选获得了398个novel lncRNA。

[novel_lncRNA.gtf](./result/novel/novel_lncRNA_end.gtf)

[novel_lncRNA_transtcript.fa](./result/novel/stringtie.fa)

#### 4. 获取基因组中已知的lncRNA并进行基因的差异表达分析

4.1 获取基因组中已知的lncRNA

> 根据基因组GRCh37.lncRNA.gtf以及stringtie_merge.gtf文件获取已知的lncRNA的注释信息(stringtie_filter.gtf)
>
> 7.filter.pl
>
> ~~~perl
> cd 6.raw.lncRNA
> #!/usr/bin/perl
> use strict;
> use warnings;
> my ($file_tie,$file_gtf) = @ARGV;
> my %gtf;
> &get_tid($file_gtf,\%gtf);
> open FILTER,">stringtie_filter.gtf" or die;
> open TIE,"$file_tie" or die;
> while (my $line = <TIE>){
>         chomp $line;
>         next if ($line =~ /^#/);
>         my $note = (split /\t/,$line)[-1];
>         $note =~ /.*transcript_id\s+\"(.*?)\";/;
>         my $tid = $1;
>         next if !exists $gtf{$tid};
>         print FILTER "$line\n";
> }
> close TIE;
> close FILTER;
> sub get_tid{
>         my ($file,$hash) = @_;
>         open OUT,">gid2tid.txt" or die;
>         open IN,"$file" or die;
>         while (my $line = <IN>){
>                 chomp $line;
>                 next if ($line =~ /^#/);
>                 my ($chr,$type,$note) = (split /\t/,$line)[0,2,8];
>                 next if ($chr !~ /chr/);
>                 next if ($type ne "transcript");
>                 $note =~ /gene_id\s+\"(.*?)\";\s+transcript_id\s+\"(.*?)\";/;
>                 my $gid = $1; my $tid = $2;
>                 $$hash{$tid} ++;
>                 print OUT "$gid\t$tid\n";
>         }
>         close IN;
>         close OUT
> }
> ~~~
>
> 共获取48056个lncRNA转录本信息
>
> [lncRNA.gtf](./result/raw/stringtie_filter.gtf)
>
> [lncRNA.fa.gz](./result/raw/stringtie.fa.gz)

4.2 利用salmon计算lncRNA的表达量

7.salmon.sh

~~~shell
cd 6.raw.lncRNA
gtf=`ls *.gtf`
$trans_path/util/gtf_genome_to_cdna_fasta.pl $gtf $ref_path/GRCh37.fa > stringtie.fa
gzip stringtie.fa

$salmon index -t stringtie.fa.gz -i stringtie.index

for i in `ls ../1.clean/*_1_val_1.fq.gz`
do
        i=${i/*\//}
        i=${i/_1_val_1.fq.gz/}
        echo $i
        $salmon quant -i stringtie.index -l A --gcBias -1 ../1.clean/$i\_1_val_1.fq.gz -2 ../1.clean/$i\_2_val_2.fq.gz -p 60 -o ./quants/$i &> ./logs/$i.log
done
~~~

利用DESeq2进行差异表达分析

Deg.R

~~~R
rm(list=ls())
workdir="~/zhaojp/lncRNA/raw_lncRNA/quants"
setwd(workdir)
library("DESeq2")
library("tximport")
library("readr")
library("ggplot2")
library('pheatmap')
genes2tx = read.table("../gid2tid.txt",sep="\t",header = TRUE)
gid = genes2tx$gid
tid = genes2tx$tid
tx2genes = data.frame(tid,gid)

samples = read.table("../samples.txt",sep="\t",header=T)

rownames(samples) = samples$run
files = file.path(".",samples$run,"quant.sf")
names(files) = samples$run
txi.tpm <- tximport(files, type = "salmon", 
                    tx2gene = tx2genes,
                    countsFromAbundance = "scaledTPM")
keep <- rowSums(txi.tpm$counts) >= 3
txi.tpm$counts <- txi.tpm$counts[keep, ]
txi.tpm$length <- txi.tpm$length[keep, ]
txi.tpm$abundance <- txi.tpm$abundance[keep, ]
dds = DESeqDataSetFromTximport(txi.tpm,
                               colData = samples,
                               design = ~ condition)
head(dds)
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "treat", "control"))
summary(res)
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
  labs(x="log2FC",y="-log10padj")+
  geom_hline(yintercept=-log10(fdr),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  ylim(0,10)+
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

pdf(file=paste("2.deg_volcano.pdf",sep=""), 
    width=8, height=8, onefile = FALSE)

print(P_volcano)
dev.off()


col = colorRampPalette(c('blue',  'white','red'))(256)
scal = "row" ####("none","row","columa")
annotation_col = data.frame(samples$condition)

rld <- rlog(dds, blind=FALSE)
eSet = assay(rld)

n.sample=ncol(eSet)
cols <- rainbow(n.sample*1.2)

pdf(file=paste("0.samples_info.pdf",sep=""), 
    width=12, height=9, onefile = FALSE)
boxplot(eSet, col = cols,main="expression value",las=2)
#p2 = hist(eSet)
dev.off()

DEG_list = rownames(deg)
DEG_list = DEG_list[order(DEG_list)]
deg_eset = eSet[match(DEG_list,rownames(eSet)),]


rownames(annotation_col) <- colnames(deg_eset)
ann_colors = list(Group=c(differentiated="green",undifferentiated="yellow"))

pdf(file=paste("1.deg_expressionHeatmap.pdf",sep=""), 
    width=9, height=12, onefile = FALSE)
p = pheatmap(deg_eset, color = col, 
         cluster_rows = T, cluster_cols=T, 
         scale = scal, show_rownames = F, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
print(p)
dev.off()

vsd <- vst(dds,blind = FALSE)
pcaData = plotPCA(vsd, intgroup = c( "condition", "type"), returnData = TRUE)
#plotPCA(vsd,intgroup = c( "condition", "type"), "batch")
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("pca.pdf",width=6, height=6)
p = ggplot(pcaData, aes(x = PC1, y = PC2, color = type, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
print(p)
dev.off()
~~~

<table><tr>
<td><img src=./result/raw/0.samples_info.pdf border=0></td>
<td><img src=./result/raw/pca.pdf border=0></td>
</tr></table>

<table><tr>
<td><img src=./result/raw/1.deg_expressionHeatmap.pdf border=0></td>
<td><img src=./result/raw/2.deg_volcano.pdf border=0></td>
</tr></table>

共获取了10个上调基因，8个下调基因：

[deg.txt](./result/raw/deg1.txt)

本次未展示GO富集结果，因为GO富集没有作出结果

