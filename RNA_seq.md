# Install

```R
# convert .sra to .fa
fastq-dump --gzip --split-3 -O ./ /public/project/RNA/airway/sra/SRR1039508.sra
ls *sra |while read id; do fastq-dump --split-3 $id;done
```

### Conda

```bash
# before install
conda search trimmomatic
# install
wget ...
bash miniconda3.sh
# set variable 
echo ". /home/<user>/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
echo ". /home/<user>/miniconda3/etc/profile.d/conda.sh" >> ~/.bash_profile
echo "conda activate" >> ~/.bashrc
export PATH="/home/<user>/miniconda3/bin:$PATH"
# update
conda update conda
# install spyder 
conda install spyder
conda install jupyter notebook
# 
conda activate
conda deactivate
#  go to check before install package https://bioconda.github.io/recipes.html
conda config --append channels https://repo.anaconda.com/pkgs/main/
# create environment
conda create -n biotools
# install package on environment
conda install -n biotools sratool
#查看环境
conda info --envs
#制作环境的完整副本
    conda create --name flowers —clone snowflakes
# 删除环境：
conda remove -n biostar -all
conda remove --name flowers --all   
conda search --full-name python
conda list | grep pandas
# remove folder
rm -rf ~/miniconda3/envs/biostar/
# active environment
source activate biotools
# to deactive an environment
source deactivate
# install package from bioconda
conda install -c bioconda fastqc 
conda install -c bioconda -y star hisat2 bowtie2
```



###Reference

```R
# form NCBI：
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz          ## hg38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/GFF/ref_GRCh37.p5_top_level.gff3.gz    ## hg19

# from ensembl：
wget ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz## hg38
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz## hg19
#################
for i in $(seq 1 22) X Y M;
do wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq/hs_ref_GRCh37.p5_chr${i}.fa.gz; 
gzip -d hs_ref_GRCh37.p5_chr${i}.fa.gz;
cat hs_ref_GRCh37.p5_chr${i}.fa.gz >>hg19.fa;
sleep 10s;
done;
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar zvfx chromFa.tar.gz
cat *.fa > hg19.fa
rm chr*.fa
# by each chromosome
for i in $(seq 1 22) X Y M;
do echo $i;
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz;
## 这里也可以用NCBI的：ftp://ftp.ncbi.nih.gov/genomes/M ... led_Chromosomes/chr前缀
done
gunzip *.gz
for i in $(seq 1 22) X Y M;
do cat chr${i}.fa >> hg19.fasta; done
rm -fr chr*.fasta
# See file:
从注释文件里找出CD28，再存到CD.gtf
grep CD28 Homo_sapiens.GRCh38.89.chr.gtf |wc > CD.gtf 
```

# QualityControl

###第一步是qc

包括使用fasqc和multiqc两个软件检查测序质量，以及使用trim_galore软件进行过滤低质量reads和去除街头。

```bash
###### fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip ./
cd FastQC
chmod 755 fastqc
/path_to_fastqc/FastQC/fastqc untreated.fq -o fastqc_out_dir/
# fastqc
fastc /home/......../SRR1039510_1.fastq.gz
fastqc -t 10 SRR1039510_1.fastq.gz (10 core)
ls *gz |xargs fastqc -t 10 ()
# combine reports
multiqc ./ 
#过滤 (切除测序接头序列和read的低质量序列)
# paired read
java -jar /path/Trimmomatic/trimmomatic-0.36.jar PE -phred33 -trimlog logfile reads_1.fq.gz reads_2.fq.gz out.read_1.fq.gz out.trim.read_1.fq.gz out.read_2.fq.gz out.trim.read_2.fq.gz ILLUMINACLIP:/path/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50
# singel read
java -jar /path/Trimmomatic/trimmomatic-0.36.jar SE -phred33 -trimlog se.logfile raw_data/untreated.fq out.untreated.fq.gz ILLUMINACLIP:/path/Trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50
### multiqc
mkdir ~/project/boy
mkdir {raw, clean, qc, align, mutation}
cd qc
mkdir {raw, clean, qc, align, mutation}
cd qc
find /../ -name *gz |grep -v '\._'|xargs fastqc -t 10 -o ./
multiqc -n row ./
```

假设质量很差，就过滤

```bash
# 去掉接头去掉质量差的数据
dir='clean'
fq1='raw/SRR1039508_1.fastq.gz'
fq2='raw/SRR1039508_2.fastq.gz'
trim_galore -q 25 --phred33 --lengh 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
# trim_galore for multiple files (paired read)
1. set variable
ls /home/jmzeng/project/airway/raw/*_1.fastq.gz >1
ls /home/jmzeng/project/airway/raw/*_2.fastq.gz >2
paste 1 2 > confid
2. write shell code
###
source activate rna
dir='/home/jmzeng/project/airway/clean'
cat $ config_file |while read id
# cat $ $1 |while read id
do
   arr=($id)
   fq1=${arr[0]}
   fq2=${arr[1]}
nohup $trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 & 
done
source deactivate
###
bash qc.sh config
```

# Align 

### Create reference index

```bash
# BWA
1 bwtsw > long genomes
bwa index -p bwa -a bwtsw hg19.fa. 
2.
nohup time ~/biosoft/bwa/bwa-0.7.15/bwa index -a bwtsw -p ~/reference/index/bwa/hg19  ~/reference/genome/hg19/hg19.fa 1>hg19.bwa_index.log 2>&1 & #
### BWA
/home/u2793/soft/bwa-0.7.17
```

### Create test file

```bash
# make test file
cd ../../test
# 如果是压缩文件
ls ../../*gz |while read id; do (zcat $id | head -1000 > $(basename $id ".gz"));done
# 如果是非压缩文件
cat fa.txt |while read id; do (zcat $id | head -1000 > $(basename $id ".gz")); done
```

### Align to reference genome (bwa)

```bash
# mapping
1.
bwa mem -t 12 -M /../hg19/bwa SHCHVD26042018_1.fastq.gz SHCHVD26042018_2.fastq.gz > bwa.sam
2. 
$ bwa bwasw genome read1.fq read2.fq > aln-pe.sam
2.
for i in $(seq 1 6) ;do (nohup ~/biosoft/bwa/bwa-0.7.15/bwa  mem -t 5 -M ~/reference/index/bwa/hg19  KPGP-00001_L${i}_R1.fq.gz KPGP-00001_L${i}_R2.fq.gz 1>KPGP-00001_L${i}.sam 2>KPGP-00001_L${i}.bwa.align.log &);done
# aligning single end reads
bwa aln -t 4 hg19bwaidx sequence1.fq.gz > sequence1.bwa
bwa samse hg19bwaidx sequence1.bwa sequence1.fq.gz> sequence1_se.sam
# aligning paired end reads
bwa aln -t 4 hg19bwaidx sequence1.fq.gz > sequence1.sai
bwa aln -t 4 hg19bwaidx sequence2.fq.gz > sequence2.sai
bwa sampe hg19bwaidx sequence1.sai sequence2.sai sequence1.fq.gz sequence2.fq.
gz > sequence12_pe.sam
# example 1 (跑gatk流程需要排序)
sample=SHCHVD26042018
bwa mem -t 12 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" /home/u2793/soft/gatk-4.0.12.0/resources/bundle/bwa_index/gatk_hg38
 HCHVD26042018_L1.fastq HCHVD26042018_L2.fastq |samtools sort -@ 12 -o $sample.bam -
# example 2
INDEX=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/bwa_index/gatk_hg38
sample=SHCHVD26042018
bwa mem -t 12 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 > $sample.sam

```

### Align in loops

```bash
## make a config file in three columns format, such as 
fq1 fq2 sample_name
#
INDEX=自己设置
cat config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
    sample=${arr[2]}
    #####################################################
    ################ Step 1 : Alignment #################
    #####################################################
    echo bwa `date`
    bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 > $sample.sam
    echo bwa `date`
    #####################################################
    ################ Step 2: Sort and Index #############
    #####################################################
    echo SortSam `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$sample.sam OUTPUT=$sample.bam 
    samtools index $sample.bam
    echo SortSam `date`
    #####################################################
    ################ Step 3: Basic Statistics ###########
    #####################################################
    echo stats `date`
    samtools flagstat $sample.bam > ${sample}.alignment.flagstat
    samtools stats  $sample.bam > ${sample}.alignment.stat
    echo plot-bamstats -p ${sample}_QC  ${sample}.alignment.stat
    echo stats `date`
    #####################################################
    ####### Step 4: multiple filtering for bam files ####
    #####################################################
    ###MarkDuplicates###
    echo MarkDuplicates `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD MarkDuplicates \
    INPUT=$sample.bam OUTPUT=${sample}_marked.bam METRICS_FILE=$sample.metrics  
    echo MarkDuplicates `date`
    ###FixMateInfo###
    echo FixMateInfo `date`
    java -Djava.io.tmpdir=$TMPDIR    -Xmx40g -jar $PICARD FixMateInformation \
    INPUT=${sample}_marked.bam OUTPUT=${sample}_marked_fixed.bam SO=coordinate  
    samtools index ${sample}_marked_fixed.bam
    echo FixMateInfo `date`
    echo ${sample}_marked_fixed.bam >>files.bamlist
    rm $sample.sam $sample.bam ${sample}_marked.bam
done 
samtools merge -@ 5  -b files.bamlist  merged.bam
samtools index merged.bam
```



### Align to reference genome with different tools 

```bash
# Hisat
1.
hisat -p 10 -x /../../genome -1 SRR1039508_val_1.fq -2 SRR1039508_val_2.fq -S tmp_hisat.sam
2.
hisat -p 10 -x $hisat_index -1 SRR1039508_val_1.fq -2 SRR1039508_val_2.fq -S tmp_hisat.sam
# Subjunc
subjunc -T 5 -i /../subjunc/hg19 -r SRR1039508_val_1.fq -R SRR1039508_val_2.fq -o tmp_subjunct.sam 
# bowtie2
bowtie2 -p 10 -x /../../genome -1 SRR1039508_val_1.fq -2 SRR1039508_val_2.fq -S tmp_bowtie2.sam
# star
```



### Generate BAM files (samtools)

```bash
#### For one file
# when SAM header present
samtools view -bS sequence1.sam > sequence1.bam 
# when no header
samtools view -bT hg19.fa sequence1.sam > sequence1.bam 
# set thread value as 12
samtools sort -O bam -@ 12 -o tmp_bwa.bam bwa.sam
# for multiple files
ls *.sam |while read if ;do (samtools sort -O bam -@ 5 -o $(basename $id ".sam").bam $id);done
# visualization
samtools view tmp_bwa.bam | less -S
# sort by coordinate to streamline data processing
samtools sort -O bam -o sequence1.sorted.bam -T temp sequence1.bam 
# a position-sorted BAM file can also be indexed
samtools index sequence1.sorted.bam 
#
samtools view BC.fg.bam
samtools view -H BC.fg.bam
```

### View bam file

```bash
samtools view bwa.bam | less -S
samtools tview bwa.bam
输入'／': 在输入chr1:143217889
samtools -reference /../hg38.fa bwa.bam
samtools flagstat bwa.bam
samtools mpileup
samtools view -H 7E5239.bam |grep iv "SQ"
```

### 最简单的找变异

```shell
# 改流程没有去除PCR重复
ref=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
time samtools mpileup -ugf $ref  example.bam | bcftools call -vmO z -o out.vcf.gz
ls example.bam |xargs -i samtools index {}
```

# GATK4

### install

```shell
wget 
mkdir -p gatk4 && cd gatk4
### 下载100G的文件才能使用
```

### bcftools

```shell
# useful links
https://www.jianshu.com/p/cc64a5663552
https://www.jianshu.com/p/49d035b121b8
cd bcftools
./configure --prefix=/home/u2793/soft/bcftools-1.9
make
make install
export PATH=/home/u2793/soft/bcftools-1.9/bin:$PATH
echo 'export PATH=/home/u2793/soft/bcftools-1.9/bin:$PATH'>> ~/.bashrc

```

### 去除PCR重复

```shell
samtools markup -r bwa.bam bwa.rm.bam
samtools markup -S bwa.bam bwa.mk.bam
```

### Markduplicate/FixMateInformation/BaseRecalibrator/ApplyBQSR/HaplotypeCaller

```shell
source activate wes
GATK=/home/u2793/soft/gatk-4.0.12.0/gatk
ref=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
snp=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/hg38/dbsnp_146.hg38.vcf.gz
indel=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

for sample in {7E5239.L1,7E5240,7E5241.L1}
do 
echo $sample  
# Elapsed time: 238.91 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
    -I $sample.bam \
    -O ${sample}_marked.bam \
    -M $sample.metrics \
    1>${sample}_log.mark 2>&1 

## Elapsed time: 442.85 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
    -I ${sample}_marked.bam \
    -O ${sample}_marked_fixed.bam \
    -SO coordinate \
    1>${sample}_log.fix 2>&1 

samtools index ${sample}_marked_fixed.bam

##  Elapsed time:411.82 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
    -R $ref  \
    -I ${sample}_marked_fixed.bam  \
    --known-sites $snp \
    --known-sites $indel \
    -O ${sample}_recal.table \
    1>${sample}_log.recal 2>&1 

# Elapsed time: 417.4 minues
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"   ApplyBQSR \
    -R $ref  \
    -I ${sample}_marked_fixed.bam  \
    -bqsr ${sample}_recal.table \
    -O ${sample}_bqsr.bam \
    1>${sample}_log.ApplyBQSR  2>&1 
    
# Elapsed time: 1,275.26 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
     -R $ref  \
     -I ${sample}_bqsr.bam \
      --dbsnp $snp \
      -O ${sample}_raw.vcf \
      1>${sample}_log.HC 2>&1  
# WES数据加-ERC GVCF参数加快运行速度。还要指定-L $bed

done

### useful links
https://broadinstitute.github.io/picard/explain-flags.html
http://www.bioinfo-scrounger.com/archives/622
```

### 分染色体载入数据

```shell
# A one-liner to split your bam file by chromosome
$ for i in {1..22} {X,Y,MT}; do echo "* "$i; mkdir -p chr/${i}; samtools view -bS <( samtools view -h in.bam $i ) > chr/${i}/out.${i}.bam; done

#针对每一条染色体都会输出一个文件夹，里面存储着数据库格式文件。
for bed in  chr{1..22} chrX chrY chrM
do
echo $bed
time $GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"   GenomicsDBImport  \
-L $bed -R $GENOME \
$(ls *raw.vcf|awk '{print "-V "$0" "}') \
--genomicsdb-workspace-path gvcfs_${bed}.db
done
# 重新整理成每条染色体的vcf文件
for bed in  chr{1..22} chrX chrY chrM
do
echo $bed 
time $GATK  --java-options "-Xmx20G -Djava.io.tmpdir=./"   GenotypeGVCFs  \
-R $GENOME  -V gendb://gvcfs_${bed}.db  -O final_${bed}.vcf
done
# 多个染色体的VCF文件合并
module load java/1.8.0_91
$GATK GatherVcfs  \
$(for i in {1..22} X Y M  ;do echo "-I final_chr$i.vcf"  ;done) \
-O merge.vcf 
```

### view VCF file

```shell
# view 查看
cat SHCHVD26042018_raw.vcf |grep -v contig|less
cat SHCHVD26042018_raw.vcf |grep -w 54586
grep -v '^##' SHCHVD26042018_filtered.vcf.recode.vcf | cut -f 1-9 | tail -10
cut -f 3 SHCHVD26042018_filtered.vcf.recode.vcf |grep -v '^ID\|#' |head -50
```

```shell
# 下载文件 CCDS.20160908.txt可以使用下面的代码：
cat CCDS.20160908.txt |perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i >exon_probe.hg38.gene.bed
```

### search variants ID

```shell
# link
https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi?view+summary=view+summary&build_id=138
```

### 统计每个样品的测序深度、覆盖度、比对率等评价指标

```shell

```

#VCF文件操作

```bash
# 有多少个变异
wc -l SHCHVD26042018.vcf
zcat /home/u2793/liuguojun/reference/gtf/gencode/gencode.v29.annotation.gtf.gz |grep -w BRCA1|head|less -SN
# 合并两个vcf 文件，合并之后记得加上第一行
for i in *.vcf ; do sed -i '1d' $i ; done
for i in *.vcf ; do cat $i >> 11.vcf; done
# 提取synnonymous
grep -v "synnonymous" variants.mm10_multianno.txt> exonic.mm10_multianno.txt
# 
awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4,$5,"Het-158"}' Het-158-indel.vcf.avinput> Het-158-indel.tsv
# 合并文件
for i in *.tsv ; do cat $i >> all.tsv; done
# 提取exnic 2
sed -n 1p tmp.hg38_multianno.txt > file1.txt
grep "exonic" tmp.hg38_multianno.txt > file2.txt
cat file1.txt file2.txt > new_exonic.txt
```

### 通过IGV可视化单个基因的BAM文件

```shell
chr17   HAVANA  gene    43044295        43170245
3.5G Jul 21 18:01 7E5240.bam
7.1G Jul 21 21:40 7E5240_bqsr.bam
4.7G Jul 21 20:28 7E5240_marked.bam
4.8G Jul 21 20:44 7E5240_marked_fixed.bam
# 
zcat /misc/home/u2793/liuguojun/reference/gtf/gencode/gencode.v29.annotation.gtf.gz |grep -w FAT4|less -SN
# then you can find the following:
chr4    HAVANA  gene    125316399       125492932
# check work count
samtools view SHCHVD26042018.bam chr4:125316399-125492932|wc
#
samtools view -h SHCHVD26042018.bam chr4:125316399-125492932|samtools sort -o SHCHVD26042018.FAT4.bam -
#
samtools view -h SHCHVD26042018_bqsr.bam chr4:125316399-125492932|samtools sort -o SHCHVD26042018_bqsr.FAT4.bam -
#
samtools view -h SHCHVD26042018_marked.bam chr4:125316399-125492932|samtools sort -o SHCHVD26042018_marked.FAT4.bam -
#
samtools view -h SHCHVD26042018_marked_fixed.bam chr4:125316399-125492932|samtools sort -o SHCHVD26042018_marked_fixed.FAT4.bam -
# 
ls *FAT4.bam|xargs -i samtools index {}
# 
mv *FAT4* tmp/
```

### 把这些bam里面的BRCA1基因的reads拿出来：

```shell
samtools  view -h  7E5240.bam chr17:43044295-43170245 |samtools sort -o  7E5240.brca1.bam -
samtools  view -h  7E5240_bqsr.bam chr17:43044295-43170245 |samtools sort -o  7E5240_bqsr.brca1.bam -
samtools  view -h  7E5240_marked.bam chr17:43044295-43170245 |samtools sort -o  7E5240_marked.brca1.bam -
samtools  view -h  7E5240_marked_fixed.bam chr17:43044295-43170245 |samtools sort -o  7E5240_marked_fixed.brca1.bam -
ls  *brca1.bam|xargs -i samtools index {}
```

### 有了这些特定基因区域的bam，就可以针对特定基因找变异

```shell
# 变异与突变不同
ref=/home/u2793/soft/gatk-4.0.12.0/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
samtools mpileup -ugf $ref   7E5240_bqsr.brca1.bam   | bcftools call -vmO z -o 7E5240_bqsr.vcf.gz
zcat 7E5240_bqsr.vcf.gz
```

### Multiqc

```bash
# 记录最多的基因是
zcat clinvar_20180429.vcf.gz|perl -alne '{/GENEINFO=(.*?):/;print $1 if $1}'|sort |uniq -c |sort -k 1,1nr >gene.clinvar.freq
head gene.clinvar.freq
# 记录的致病情况最多的基因是
zcat clinvar_20180429.vcf.gz|perl -alne '{/GENEINFO=(.*?):/;$g=$1;/CLNSIG=(.*?);/;print "$g\t$1" if $1}'| sort |uniq -c |sort -k 1,1nr >gene_sign.clinvar.freq
head gene_sign.clinvar.freq
# 看看是否有被专家审核
zcat clinvar_20180429.vcf.gz|perl -alne '{/GENEINFO=(.*?):/;$g=$1;/CLNSIG=(.*?);/;$t=$1;/CLNREVSTAT=(.*?);/;print "$g\t$t\t$1" if $1}'| sort |uniq -c |sort -k 1,1nr >gene_sign_review.clinvar.freq
[jianmingzeng@jade anno]$ head gene_sign_review.clinvar.freq
```

### cut

```shell
### cut file
grep "#CHROM" genotype_id.vcf |cut -f 10-50 |tr '\t' '\n' > sample_id.txt
```

#Filtering  variants in vcf file with VCFtools

### install

```shell
# link
群体遗传专题：关于SNP的过滤（1）：如何使用vcftools进行SNP过滤
http://www.ddocent.com/filtering/
# install
git clone https://github.com/vcftools/vcftools.git
cd vcftools
mkdir /home/u2793/soft/vcftools -p 
./autogen.sh
./configure --prefix=/home/u2793/soft/vcftools
make
make install
echo 'export PATH=/home/u2793/soft/vcftools/src/cpp:$PATH'>> ~/.bashrc
source ~./bashrc
```

### 过滤标准

```shell
# For SNPs:
QD<2.0
MQ<40.0
FS>60.0
SOR>3.0
MQRankSum<-12.5
ReadPosRankSum<-8.0
# For indels:
QD<2.0
ReadPosRankSum<-20.0
InbreedingCoeff<-0.8
FS>200.0
SOR>10.0
# for WGS
# 第8列的DP4进行过滤：vcf文件中的 DP4 一列信息非常重要，提供了4个数据：1：比对结果和正链一致的reads数、2：比对结果和负链一致的reads数、3：比对结果在正链的variant上的reads数、4：比对结果在负链的variant上的reads数。可以设定 （value3 + value4）大于某一阈值，才算是variant。比如：
perl -ne 'print $_ if /DP4=(\d+),(\d+),(\d+),(\d+)/ && ($3+$4)>=10 && ($3+$4)/($1+$2+$3+$4)>=0.8' smaple_raw.vcf > sample_filter.vcf

```

### usage

```shell
###已存在于数据库中的变异写入 *dropped文件，在数据库中不存在的变异信息将会被写入到*filtered
# use
vcftools --vcf SHCHVD26042018_raw.vcf --minDP 10 --minQ 30 --recode --recode-INFO-all --out SHCHVD26042018_filtered
vcftools --vcf SHCHVD26042018_raw.vcf --minDP 10 --minQ 30 --recode --out SHCHVD26042018_filtered
# 删除任何indel位点
vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only
# 没有任何具有过滤器标记的位点，然后使用gzip压缩它
vcftools --gzvcf input_file.vcf.gz --remove-filtered-all --recode --stdout | gzip -c > output_PASS_only.vcf.gz
#从输入vcf文件输出新的vcf文件，该文件删除任何indel位点
vcftools --gzvcf input_file.vcf.gz --freq --chr 1 --out chr1_analysis
#
vcftools --gzvcf raw.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3
# 测可能产生的错误率(小于0.05有意义)
curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ErrorCount.sh
chmod +x ErrorCount.sh 
./ErrorCount.sh raw.g5mac3dp3.recode.vcf 
#
vcftools --vcf raw.g5mac3dplm.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out DP3g95maf05 --min-meanDP 20
#
vcftools --vcf genotype_id.vcf -maf 0.05 --min-alleles 2 --max-alleles 2
# link
https://www.jianshu.com/p/e05ff3cace56
https://mp.weixin.qq.com/s?src=11&timestamp=1553615209&ver=1508&signature=PFuITkp4ZBXDX7ORhhBPg-Ctf3dS80z0F1eWsMEFnU4JJLcqDuB*cTp4xWZBiCHkbZ6pqZk9PnM7d48qG0vUsYaiY4Mm*S3KMk*pm6-DqTx27BREtnZAHX9uwSJCE-1T&new=1
https://mp.weixin.qq.com/s?src=11&timestamp=1553615209&ver=1508&signature=PFuITkp4ZBXDX7ORhhBPg-Ctf3dS80z0F1eWsMEFnU60EuRnRYrIo0pqQX3U85Oi6jdTWCbQVRUWtdhQpcZkUsq53EALPOCYaIkmjRxXIshqKAsanXfjidg0ShtvMMzs&new=1
https://www.jianshu.com/p/52b2dcb601d2
```

#ANNOVAR

### 下载数据库

```shell
# usefull link 
http://blog.applymed.cn/index.php?k=NGS&id=3_5
http://www.bio-info-trainee.com/3838.html
http://annovar.openbioinformatics.org/en/latest/user-guide/download/
# download one database
/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 refGene /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 knownGene /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 avsnp150 /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 intervar_20180118 /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 cosmic70 /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 exac03 /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 1000g2015aug /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 gnomad_genome /home/u2793/soft/annovar/humandb/

/home/u2793/soft/annovar/annotate_variation.pl -downdb -webfrom annovar -buildver hg38 dbnsfp35c /home/u2793/soft/annovar/humandb/

# download database in loops
DATABASE_LIST="
refGene
avsift
ljb26_all
cosmic68wgs
cosmic70
esp6500siv2_all
1000g2010nov
1000g2014oct
snp131
snp138
snp131NonFlagged
snp138NonFlagged
clinvar_20150629
"
for DATABASE in $DATABASE_LIST
do
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar $DATABASE humandb/
done
# cytoBand AND genomicsSuperDups need to be downloaded seperately
./annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
./annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
```



### 根据质量值过滤

```shell
#GATK variantfiltration
#vcftools
# link
https://mp.weixin.qq.com/s/W8Vfv1WmW6M7U0tIcPtlng
```



### 基于基因的注释

```shell
#
~/soft/annovar/convert2annovar.pl -format vcf4old SHCHVD26042018_raw.vcf >SHCHVD26042018.annovar
#
~/soft/annovar/annotate_variation.pl \
-buildver hg38  \
--outfile SHCHVD26042018.anno \
SHCHVD26042018.annovar \
~/soft/annovar/humandb/
# search gene
cat SHCHVD26042018.anno.exonic_variant_function |grep -w F7
# check variants on exon 
wc -l SHCHVD26042018.anno.exonic_variant_function
less SHCHVD26042018.anno.exonic_variant_function
# 处理掉那些无义突变
sed 's/^NA/unknown/' var_annovar_maf > var_annovar_maf2
grep -v "^NA" var_annovar_maf | grep -v -P "\tUNKNOWN\t"> var_annovar_maf2
# 兴趣的snp位点（dbsnp rsID）且只想对这些位点进行注释
cat example/snplist.txt 
rs74487784
rs41534544
rs4308095
rs12345678
convert2annovar.pl -format rsid example/snplist.txt -dbsnpfile humandb/hg19_snp138.txt > snnplist.avinput
```

### Region-based annotation

```shell

```

### 基于数据库的注释

```shell
#目标变异筛选(基于变异频率)
结合以上数据库，通过特定的阈值筛选，我们可以过滤很多无效变异。例如，可以过滤千人基因组数据库中频率大于0.01变异位点，以得到真正可能致病的罕见突变（rare）。也可以联合多个数据库对突变频率进行过滤，或者同时参考dbSNP中记录的SNP信息，初步判断数据库中不存在的变异为新发现变异，以增加研究价值。
# 重点是检查那些后缀为 dropped 的文件
# vcf 先转成annovar file
perl convert2annovar.pl -format vcf4 cohort.vcf > cohort.avinput
# vcf 先转成annovar file (基于多样本)
perl /Users/chengkai/Documents/05_software/annovar/annovar/convert2annovar.pl -format vcf4 -allsample -withfreq  cohort.filter.vcf > cohort.filter.avinput
# 使用1000Genomes数据库进行频率注释
annotate_variation.pl -filter -dbtype 1000g2012apr_eur -buildver hg19 -out ex1 example/ex1.avinput humandb/
# 手动制作文件
cut -f 5-7,12,13,1,16 human_brca_all_mutect2.maf |cut -f 2-7  > 1
cut -f 5-7,12,13,1,16 human_brca_all_mutect2.maf |cut -f 1 > 2
paste 1 2 > for_annovar.input ### 共 13027 位点
# 筛选对蛋白功能影响的snvs和indel(sift, polyphen-2, provean,dbSNP132，1000GenomeProject，HapMap)
~/soft/annovar/annotate_variation.pl -filter -dbtype dbnsfp35c -buildver hg38 -out ex1 SHCHVD26042018.anno.exonic_variant_function ~/soft/annovar/humandb annotate_variation.pl -filter -out ex1 -build hg19 -dbtype snp138 example/ex1.avinput humandb/
# 过滤
-score_threshold
--normscore_threshold
-maf 0.05 -reverse
-maf 0.05
# 筛选与疾病，临床，免疫相关的(clinvar,HGMD)
# 筛除/排除/过滤常见的变异(用ExAC，dbSNP，千人基因组计划数据库)
#filter rare or unreported variants (in 1000G/dbSNP) or predicted deleterious variants
# 注释多个数据库
# csv output file 
/home/u2793/soft/annovar/table_annovar.pl SHCHVD26042018_raw.avinput /home/u2793/soft/annovar/humandb/ -buildver hg38 -out cohort -remove -protocol refGene,cytoBand,dbnsfp35c,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp150,clinvar_20190305,intervar_20180118,gnomad_genome -operation g,r,f,f,f,f,f,f,f,f,f -nastring . -csvout
# .vcf output file
/home/u2793/soft/annovar/table_annovar.pl SHCHVD26042018_raw.avinput /home/u2793/soft/annovar/humandb/ -buildver hg38 -out cohort -remove -protocol refGene,cytoBand,dbnsfp35c,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp150,clinvar_20190305,intervar_20180118,gnomad_genome -operation g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
# .txt output file (default is .txt)
bin=/home/u2793/soft/annovar/
db=/home/u2793/soft/annovar/humandb/ 
$bin/table_annovar.pl SHCHVD26042018_raw.avinput $db -buildver hg38 -out  tmp \
-protocol refGene,cytoBand,dbnsfp35c,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp150,clinvar_20190305,intervar_20180118,gnomad_genome \
-operation g,r,f,f,f,f,f,f,f,f,f -nastring NA
# 
/home/u2793/soft/annovar/table_annovar.pl SHCHVD26042018_raw.avinput /home/u2793/soft/annovar/humandb/ -buildver hg38 -out cohort -remove -protocol refGene,cytoBand,dbnsfp35c,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp150,clinvar_20190305,intervar_20180118,gnomad_genome -operation g,r,f,f,f,f,f,f,f,f,f -nastring . -csvout -polish -xref gene_fullxref.txt

# example 2
perl  cohort.filter.avinput /Users/chengkai/Documents/05_software/annovar/annovar/humandb/ -buildver hg19 -out cohort -remove -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,exac03,avsnp147,dbscsnv11,clinvar_20180603,intervar_20170202,icgc21,dbnsfp33a,gnomad_exome -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -polish

### 使用table_annovar.pl注释多个数据库
#-bulidver hg19 表示使用的参考基因组版本
#-out myanno 表示输出文件前缀，亦可用 --outfile 直接指定文件名
#-remove 表示删除中间文件
#-protocol 表示使用的数据库，其数据库顺序要与后面的operation注释方式对应上
#-operation 表示对应数据库的注释类型（g代表gene-based、r代表region-based、f代表filter-based，gx means gene-based with cross-reference annotation (from -xref argument)）
#-nasting . 点号代替缺省值
#-csvout 表示输出为csv格式
# -polish, It polish the amino acid predictions for indels, and can be very useful to know the exact change to protein sequence, rather than just a simple "frameshift" annotation.
# 除了1000g之外，其他库名都要和文件名一致
perl table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
# link
https://www.jianshu.com/p/7607c894eaae
```

### 去除高频率变异

```perl
# link
https://blog.csdn.net/weixin_33895604/article/details/87201102
https://www.jianshu.com/p/64b15ab187ab?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation
perl -alne '{print if $F[1]>0.05}' tmp.hg38_ALL.sites.2015_08_dropped > filter_by_1000g.pos
cat filter_by_1000g.pos  ping_organoids_all.maf   |perl -alne '{if(/^1000/){$h{"$F[2]\t$F[3]"}=1}else{print unless exists $h{"$F[4]\t$F[5]"}}}' > filter_by_1000g.maf
perl -alne '{print if (split(",",$F[1]))[0]>0.05}' tmp.hg38_gnomad_genome_dropped > filter_by_gnomad.pos
perl -alne '{print if (split(",",$F[1]))[0]>0.05}' tmp.hg38_exac03_dropped  > filter_by_exac03.pos
cat filter_by_exac03.pos filter_by_1000g.maf |perl -alne '{if(/^exac03/){$h{"$F[2]\t$F[3]"}=1}else{print unless exists $h{"$F[4]\t$F[5]"}}}' >  filter_by_1000g_exac03.maf
cat filter_by_gnomad.pos   filter_by_1000g_exac03.maf|perl -alne '{if(/^gnomad_genome/){$h{"$F[2]\t$F[3]"}=1}else{print unless exists $h{"$F[4]\t$F[5]"}}}' >  filter_by_1000g_exac03_gnomeAD.maf
wc -l *.maf
# link
http://www.biotrainee.com/thread-1376-1-1.html
# ANNOVAR人类各个数据库变异注释结果表格说明
http://www.omicsclass.com/article/464
https://brb.nci.nih.gov/seqtools/colexpanno.html
```

### 倒入R

```R
####Rstudio：导入到maftools###
library(maftools)
var.annovar.maf = annovarToMaf(annovar = "variants.mm10_multianno.txt",Center = 'CSI-NUS', refBuild = 'mm10',tsbCol = 'Tumor_Sample_Barcode', table = 'refGene')
write.table(x=var.annovar.maf,file="var_annovar_maf",quote= F,sep="\t",row.names=F)
var_maf = read.maf(maf="var_annovar_maf")
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median',showBarcodes = T)
oncoplot(maf = var_maf, top = 10, fontSize = 12)
```



### WGS links

```shell
https://blog.csdn.net/weixin_33893473/article/details/87561225
```

### search SNP (rs)

```shell
https://www.ncbi.nlm.nih.gov/snp/rs788908#seq_hash
```

