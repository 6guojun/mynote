 

# Bash

### Basic

```shell
# help
man grep
# wc
wc -l: 行数
wc -w: 字数
wc -c: 字符数
# tar
用于多个文件或目录进行打包，但不压缩; 也用于解包
-c 创建新文档
-x 解包
-v 显示正在处理的文件名
-f 取代默认的文件名
# gzip, 用于文件进行压缩和解压缩命令，文件扩展名为.gz结尾
# gunzip, 用于对gzip压缩文档进行解压缩
# zip/unzip, 压缩解压缩.zip文件
打包： tar -cvf xxx.tar *
压缩： gzip xxx.tar
打包： gzip -d xxx.tar.gz
压缩： tar -xvf xxx.tar
打包并压缩：tar -zcvf xxx.tar.gz *
解压：tar -zxvf xxx.tar.gz
打包并压缩：tar -jcvf xxx.tar.bz2 *
解压：tar -jxvf xxx.tar.bz2
# 查看ip地址
ifconfig：
# 查看某列，column -t 是为了好观察，不要把用它处理然后把结果传入文本
grep -v "^#" Homo_sapiens.GRCh37.75.gtf | head -n 10 | cut -f 3-5 | column -t
# -c显示计数
sort test.letter | uniq -c
# -d选项只输出重复行
uniq -d test.letter
```



```shell
#
date; who
# to check cpu model
sysctl -n machdep.cpu.brand_string
#History:
history
# Display current path/direction:
pwd
dirname test.txt
#space
df -h
# top big folders
du -a | sort -n -r | head -n 20
du -a /home | sort -n -r | head -n 5
du -sh /Users/liu/Desktop/*
# Search direction
which fastqc
#To change directory：
cd desktop
cd .., cd ~, cd -, cd /
cd ../../bowtie2 
go to hard disk:
cd /Volumes
# lists your files in 'long format'
ll 
ls -lh
ls -al
ls -lh |wc
ls raw/
ls raw/*
ls
ls -al: 
ls -1 
ls -lh
ls -lh *gz
ls */*/
ls ../clean/*gz
ls *gz |while read id;do (zcat $fwelkjfewk|head -1000);done
# To Get help:
man mkdir
mkdir -h

type ls查看是不是内置命令
# check server version
lsb_release -a
check IP configuration:
      ifconfig
To turn into next line:
     \return
To check network (internet):
      ping example.com
      ping -c 1 youtube.com
      ping -c 1 example.com | cut -d = -f 4
      echo `ls -al`:优先执行 ls -al,再打印。
# To list files/show the files
#Count:
 wc 
 wc -l
# See file size
du -h
du -h SHCHVD26042018_1.fastq
# 文件夹操作：
See folder
    echo $PATH | tr ‘:’ ‘\n’
    echo folder{1..5}/folder{1..5} | xargs -n 1
# Directory/folder
      mkdir tmp
      mkdir -p 1/2/3
      mkdir {raw,clean,qc,align,mutation}
      mkdir {1..10}
      mkdir {1..10}/{1..10}
      mkdir samtools && cd samtools
      mkdir -pv Folder_{1..5}/{1..5}
# proseccer
ps
# kill the proseccer
kill -9 id
#stop the process
ps -ef | grep prefe |awk '{print $2}' | while read id | do kill $id ; done
# 添加可执行权限并且执行该脚本
chmod u+x test13
./test13
```

### Variable

```shell
# see path
echo $PATH
vi ~/.bashrc
#Environment 
#进入文件夹
cd Users/liu//biosoft/bowtie2/bowtie2-2.3.4.3-macos-x86_64
./bowtie2
#设置变量
bowtie2=/Users/liu//biosoft/bowtie2/bowtie2-2.3.4.3-macos-x86_64/bowtie2
$bowtie2
#设置alias(临时)
alias bowtie2=‘/Users/liu//biosoft/bowtie2/bowtie2-2.3.4.3-macos-x86_64/bowtie2’
#设置环境变量(临时)
$export PATH=“$PATH:/Users/liu/biosoft/bowtie2/bowtie2-2.3.4.3-macos-x86_64/”
#在bashrc设置环境变量(永久)
cat >> ~./bashrc 
source ~./bashrc (该文件source后生效)
# xargs
ls *.sam | xargs -i head {}
# 映射
ln -s /usr/local/mysql/bin/ mysql
```

### Unzip

```shell
# zip file
unzip file.zip
tar zxvf 文件名
#.tar.bz2
tar -jxvf samtools-1.8.tar.bz2
#.gz file       
gunzip *.gz
gzip -d SHCHVD26042018_1.fastq.gz
gunzip file.fastq.gz
time gunzip file.fastq.gz
# unzip in backgroud
 & (e.g. tar -czf home.tar.gz &)
ls *zip|while read id;do unzip $id;done
# see zip file
zcat
zless
```

### echo

```shell
echo This is a test

echo "let's go"
```

### expression

```shell
# 它只支持整数运算！
var1=$[1+5];var2=$[3];re=$[var1*var2];echo $re
# 
bc
12*5.2
3.156 * (3 + 5)
quit
#
echo "scale=3; 3.44/5" | bc
#
var1=10.46
var2=43.67
var3=33.2
var4=71
var5=$(bc <<EOF
scale=4
a1 = ( $var1 * $var2)
b1 = ( $var3 * $var4)
a1 + b1
EOF
)
echo $var5
```



# Text precess

```R
# 将文件的某一列切出来
cut -f 1,3,5 filename 获取文本的第一、三、五列
cut -f 1-5 filename 获取文本的第一至五列
文件名分割
ls *gz |cut -d"_" -f 1 |sort -u
# 先按/分段后取第一个字段是a
echo "a/b/c" |cut -d '/' -f 1
cut -d '/' -f 1
# 显示行号 
nl
# 去重复
cat filename | cut -f 2 | sort | uniq | nl
```

### Read and write

```shell
Write file
TXT:
     less -S mRNA.txt
     
     head -n1 example_gene_annotations.txt
Number of rows:
     cat BLCAmRNA.csv | wc -l
     wc -l filename.csv
     wv -l *


Create File:
     touch gene.txt
     touch {1..5}/{1..5}/tmp.txt
     vim test.txt
     cat >test.txt (可以直接输入)
```

### Copy, remove, delete

```shell
# Copy
cp -r
cp -i /home/levan/kdenlive/untitelds.mpg /media/sda3/SkyDrive/
# copy
cp file.doc newfile.doc
# rename
rm file.doc newfile.doc
# move
mv *.jpg /Users/liu/Documents/01Wallpaper/
# remove
rm -rf finename
rm gene*
# remove directory 
rm -rf myfile # empty folder
rm -r myfile # non-empty folder
# Empty the trash:
rm -rf ~/.Trash/* 
```

```bash
head -1000 test.bed							
tail -10 test.bed
more test.bed
# see tail
less -S TAIR10_GFF3_genes.gff |head
less -SN test.bed
zless -SN test.fastq.gz
# no to see N
grep -v “NNN”
# See rows
     wc -l *
     cat test.bed
# See Columns：
# 截取列1到3列
cat NC.gff | cut -f 1,2,3 | head -5
     cut -f 2 Homo_sapiens.GRCh38.89.chr.gtf 
     cut -f 2 Homo_sapiens.GRCh38.89.chr.gtf |sort -u
     cut -f 1-3 test.bed
     查看某列有多少基因
     cut -f 9 Homo_sapiens.GRCh38.89.chr.gtf |cut -d";" -f 6|sort -u |wc
# write file:
      Cat:
      cat gene.txt
      cat > tmp.txt
      cat >> ~./bashrc
      cat -n test.bed(添加列的顺序)
      nano: nano.bash_profile
三驾马车：
     grep -n TP53 test.bed
# Sort
按第二列反向排序：
sort -k2,2nr test | cut -f 1-3
# -k1,1V，自然 chr11不会在chr2前面
sort -k1,1V -k2,2n test.bed
Configure terminal theme:
    PS1="\[\e[32;1m\]\u \[\e[33;1m\]\t \[\e[35;1m\]\w \n\[\e[0;40m\]$"
    export CLICOLOR=1
 export LSCOLORS=ExFxBxDxCxegedabagacad
 alias ls='ls -GFh'
Run bash:
       bash -x test.sh

Open anything:
     open -a Pages.app name.doc
     open -a Preview 编程大全.第3版.pdf
Edit File:
     nano 1064genes.csv
字符串操作
~/str/:~ 表示模式开始。/ /中是模式.    
~ /str1|str2/:是匹配str1或者str2.
!~/str/: 除了str都匹配  ，或是’!／str／’
grep -v：排除
awk '{print FILENAME"\t"$0}' * |grep -v EnsEMBL_Gene_ID >tmp.txt
#打印 “文件名+tab+全部内容”，注意排除EnsEMBL那行，输入到tmp.txt 
在左边插入新的两列：
awk '{$0=(NR==1? "NAME\tnote":FILENAME"\t") "\t" $0}7' file 
input:
head1 head2 head3
value1 value2 value3
Output:
NAME    note    head1 head2 head3
file            value1 value2 value3
#多动作
Download
     wget -c http://www.biotrainee.com/jmzeng/igv/test.bed
显示中文：
iconv -c -f GB2312 -t UTF-8 矩阵转置.cpp >> 矩阵转置2.cpp
control+d　　　　　中断a.out运行
nano 　　　　　　编写脚本语言　　ctrl+o存储
nano ....sh　　　　打开

echo "...$i..."　　　输出语句

To search in text:
     grep -i break-in auth.log
     grep stage gene.txt | awk {‘print $2’}
     grep -no H3K4me3 test.bed
     grep -n --color  H3K4me3 test.bed
awk:    ’{}’
$0表示行,$1……$n表示列
$0:所有内容
$1~$n:当前记录的第n个字段，字段间由FS分隔
FS：输入字段分隔符 默认是空格或Tab.  
-F:输入字段分隔符 默认是空格或Tab. 
NF：当前记录中的字段个数，就是有多少列
NR：Number of Records，行号，从1开始
FNR：当前记录数，与NR不同的是，这个值会是各个文件自己的行号
RS：输入的记录分隔符， 默认为换行符
OFS：输出字段分隔符， 默认也是空格
ORS：输出的记录分隔符，默认为换行符
FILENAME：当前输入文件的名字
字符串匹配：
~/str/:~ 表示模式开始。/ /中是模式.    
~ /str1|str2/:是匹配str1或者str2.
!~/str/: 除了str都匹配  ，或是’!／str／’
grep -v：排除
实际例子：
多个文件合并成一列。
awk '{print FILENAME"\t"$0}' * |grep -v EnsEMBL_Gene_ID >tmp.txt
#打印 “文件名+tab+全部内容”，注意排除EnsEMBL那行，输入到tmp.txt 
在左边插入新的两列：
awk '{$0=(NR==1? "NAME\tnote":FILENAME"\t") "\t" $0}7' file 
input:
head1 head2 head3
value1 value2 value3
Output:
NAME    note    head1 head2 head3
file            value1 value2 value3
for g in *.gz; do gunzip $g; done
for f in *.gz ; do gunzip -c "$f" > /home/$USER/"${f%.*}" ; done
/Users/liu/Documents/01Data/01TCGA/01KIRC/gdc_download_from_gdc/00_data_read_in_one_file/
Calculate:
       bc
       echo "$((5 * 5))”
       function calc { bc -l <<< ${@//[xX]/*}; };
              calc 5+1
              calc 3/1
              calc 3x2
      echo "3 * 2.19" | bc -l
      expr 1 \* 3
     expo 1 / 5
String Manipulation

Text

Set symbol between two  columns:
      awk -F':' 'BEGIN{OFS="=";} {print $3,$4;}’ gene.txt
Divide context:
       awk ‘/^sample/{n++}{print > “output.txt”}’ gene.txt
Json:
     curl https://jsonplaceholder.typicode.com/posts
     curl -i https://jsonplaceholder.typicode.com/posts/2
     curl -o test.txt https://jsonplaceholder.typicode.com/posts
    curl -O https://jsonplaceholder.typicode.com/posts
     
       loops
for i in 1 2 3 4 5
do 
    echo $i 
done

for i in {1..100}
do 
    echo $i 
done

for i in {1..100..2}
do 
    echo $i 
done
# 拆分
split -2 README #将README文件每2行分割成一个文件
-b 字节数： 按照字节数切分
split -b 200 infile
```

### join

```shell
# 两个文件都以第一列为标准
join -1 1 -2 1 example_sorted.bed example_length_alt.txt
# -a选项指定哪一个文件可以不遵循配对
join -1 1 -2 1 -a 1 example_sorted.bed example_length_alt.txt
```

# sed

### show

```shell
# link 
https://www.cnblogs.com/ctaixw/p/5860221.html
# 例如,输出偶数行， 
sed ‐n '2~2 p' test.txt
输出奇数行内容。
sed ‐n '1~2 p' test.txt
# 打印第一行到第三行
sed -n '1/3p' data
# 打印匹配"second"的行
sed -n '/second/p' data
# 打印包含first的行到第四行
sed -n '/first/, 4p'
# 打印匹配hello, world的所有行
sed ‐n '/hello/, /world/ p' test.txt
###判断文件时制表符还是空格分割的；
sed -n 1 sampleID.txt 
# 将制表符修改为逗号分隔符：
cat sampleID.txt|tr  "\t" "," > sampleID_form.txt　　
# first line, second
sed -n 1p tmp.hg38_multianno.txt
sed -n 2p tmp.hg38_multianno.txt
# 输出关键字ruby所在行的内容
sed -n '/ruby/p' ab.txt
# 输出关键字$所在行的内容，使用反斜线\屏蔽特殊含义
sed -n '/\$/p' ab.txt
```



### replace

```shell
sed 's/ruby/bird/g' ab.txt   # 把全部的ruby替换为bird
sed 's/chrom/chr/' a.txt # replace,将chrom替换成chr
sed 's/ruby//g' ab.txt   # 把全部的ruby替换为空，即删除ruby字符串
```



```shell
# 所有的命令都会一个叫做在模式空间（pattern buffer）的缓冲区进行。因此不会改变原始输入文件的内容
# link
http://dongweiming.github.io/sed_and_awk/#/
sed:编辑
sed - - help
# Number of columns:
head -1 BLCAmRNA.csv |sed 's/[^,]//g' |wc -c
# Join two file (join要求文件必须排序,所以使用sort命令,而我们需要按照两个文件的第一列合并,所以是sort -k 1,1)
join <(sed '/^#/d' a.txt|sort -k 1,1) <(sed '/^#/d' b.txt|sort -k 1,1)
# m,+n表示从m行开始向下n行， m~n表示从m行开始的每n行。

sed '1d' ab.txt         # 输出删除第一行后的文件内容
sed '$d' ab.txt         # 输出删除最后一行后的文件内容
sed '1, 2d' ab.txt      # 输出删除第一行到第二行后的文件内容
sed '2, $d' ab.txt      # 输出删除第2行到最后1行后的文件内容
# 搜索search
# 在第一行后增加字符串"drink tea"
sed '1a drink tea' ab.txt     
# 在第一行到第三行后增加字符串"drink tea"
sed '1,3a drink tea' ab.txt  
# 在第一行后增加两行，换行使用\n，可多次使用\n添加多行
sed '1a drink tea\nor coffee' ab.txt 
# 行替代
sed '1c Hi' ab.txt    # 把ab.txt的第一行替换为Hi
sed '1,2c Hi' ab.txt  # 把ab.txt的第一行到第二行替换为Hi
# first two lines
sed -n 1,2p tmp.hg38_multianno.txt
awk 'FNR <= 2' tmp.hg38_multianno.txt
perl -ne'1..2 and print' /etc/passwd
# 对每行匹配到的第一个字符串进行替换
sed -i 's/原字符串/新字符串/' ab.txt 
# 对全局匹配上的所有字符串进行替换
sed -i 's/原字符串/新字符串/g' ab.txt 
# 删除所有匹配到字符串的行
sed -i '/匹配字符串/d'  ab.txt  
# 特定字符串的行后插入新行
sed -i '/特定字符串/a 新行字符串' ab.txt 
# 特定字符串的行前插入新行
sed -i '/特定字符串/i 新行字符串' ab.txt
# 把匹配行中的某个字符串替换为目标字符串
sed -i '/匹配字符串/s/源字符串/目标字符串/g' ab.txt
# 在文件ab.txt中的末行之后，添加bye
sed -i '$a bye' ab.txt   
# 对于文件第3行，把匹配上的所有字符串进行替换
sed -i '3s/原字符串/新字符串/g' ab.txt 


```

# grep

```shell
# 注释人间中查找某个文件
grep "AT5G25475" TAIR10_GFF3_genes.gff  | head -5
# 希望搜索AT5G25475但是不包含CDS
grep "AT5G25475" TAIR10_GFF3_genes.gff  | grep -v 'CDS' | head -5
# 基因组文件中，查找某一段特定序列并查看上下文(-A n 显示后n行，-B n 显示前n行，-C n显示前后n行)
grep -A 2 'TTATTGTTGTTAAGAAAAAAGG' TAIR10_chr_all.fa
# 某一段序列出现的次数
grep  'TTATTGTTGTTAAGA' TAIR10_chr_all.fa  | wc -l
grep  'TTATTGTTGTTAAGA' TAIR10_chr_all.fa  -c
# 只返回匹配到的内容
grep  'TTATTGTTGTTAAGA' TAIR10_chr_all.fa  -o
# 
grep 文本搜索查找与匹配
无参数选项： 直接加查找内容:  grep "查找内容" filename
grep '查找内容' filename --color=always：表示突出显示查找内容
参数选项：
grep -i pattern files ：不区分大小写地搜索。默认情况区分大小写，
grep -l pattern files ：只列出匹配的文件名，
grep -L pattern files ：列出不匹配的文件名，
# grep -w，完全匹配，只匹配整个单词，而不是字符串的一部分(如匹配’magic’，而不是’magical’)，

grep -C number pattern files ：匹配的上下文分别显示[number]行，
grep pattern1 | pattern2 files ：显示匹配 pattern1 或 pattern2 的行，
grep pattern1 files | grep pattern2 ：显示既匹配 pattern1 又匹配 pattern2 的行。
grep -n pattern files 即可显示行号信息
# grep -c pattern files 即可查找总行数，即某一匹配内容出现的个数
grep -c "[^ \\n\\t]" new_exonic.txt
grep -v pattern files 查找不包含匹配项的行
# grep -o pattern files 只返回匹配到的内容，-E指定支持扩展表达式
grep -E -o 'gene_id "(\w+)"' Homo_sapiens.GRCh37.75.gtf | cut -f2 -d" "| sed 's/"//g' | sort | uniq | head -n 10

# 如： 
cat filename | grep -v "不匹配内容' | cut -f 2,5,9 | wc -l
如果想查找一个以AT5G254开头以1结尾的基因，要用到强大的正则表达式。grep 'AT5G254.*5$' TAIR10_GFF3_genes.gff
```

# 

# awk

### print

```shell
# $0 表示所有列，$1，$2`...等等表示对应的列
# print 
awk '{ print $2 "\t" $3}' example.bed
awk '{ print $2 $3}' example.bed
awk '{ print $2 , $3}' example.bed
# 第一列的第四第五字符
awk '{print $1}' new_exonic.txt | cut -c4,5
# 打印时，第一列前面加“chr”字符
awk '{ print $1}' example.bed | cut -c4 | awk '{print "chr"$1}'
# 加“Het-158”列
awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4,$5,"Het-158"}' Het-158-indel.vcf.avinput> Het-158-indel.tsv
# 这里 ~ 符号用来匹配正则表达式 
awk '$1 ~/chr1/ && $3 - $2 > 10' example.bed
# 查看.gtf 
head -n1000 Homo_sapiens.GRCh37.75.gtf | awk '!/^#/{ print $1 "\t" $4-1 "\t" $5} ' | head -n 3
# FS  输入列分割符
OFS  输出列分隔符
RS  输入行分割符
ORS  输出行分割符
NF  记录数（当前处理的）
NR  行数（当前处理）
FNR  行数（当前文件）
```



```shell

# 先把之前的tab分隔文件弄成逗号分隔文件，再查看
grep -v "^#" Homo_sapiens.GRCh37.75.gtf | head -n 10 | cut -f 3-5 | awk '{FS="\t";OFS=",";}{print $1,$2,$3}'
# seperate columns using |, OFS reffers to output file seperate
awk 'BEGIN{OFS="|"} {print $1,$2,$3}' file.txt
# show THREE columns 5 rows, $1, $2, $3指的就是第1-3列的数据
cat NC.gff | awk '{print $1,$2,$3}' | head -5
# Awk
cut -f 1-3 test.bed | awk ‘{print $1”:”$2”,”$3}’
# Number of rows and columns for CSV:
awk -F, 'END {printf "Number of Rows : %s\nNumber of Columns = %s\n", NR, NF}' BLCAmRNA.csv
# 看第二列
awk '{print $2}' abc.txt
# 第5列的数值减去第4列的数值后+1
cat NC.gff | awk ' { print $3, $5-$4 + 1 } ' | head -5
# 计算所有基因的累积长度。
cat NC.gff | awk '$3 =="gene" {len=$5-$4 + 1; size += len; print "Size:", size }'
# 显示所有的列
awk '{print $0}' TAIR10_GFF3_genes.gff | head -2
# 显示1，4，5列
awk '{print $1,$4,$5}' TAIR10_GFF3_genes.gff | head -2
# 依次显示4，5，1 列， 并且用“，”隔开
awk '{print $4","$5","$1}' TAIR10_GFF3_genes.gff | head -1
# 找到长度大于10kb且在一号染色体的注释内容
awk '$5 - $4 > 10000 && $1 ~ /Chr1/' TAIR10_GFF3_genes.gff  | head -5
# awk还有两个特殊模式BEGIN,END，顾名思义就是在操作开始或/和结束后才执行的操作。
awk 'BEGIN  {s = 0;line = 0 } ;
$5 - $4 > 10000 && $1 ~ /Chr1/ { s += ( $5 - $4 );
line += 1}; 
END {print "mean=" s/line}' TAIR10_GFF3_genes.gff  | head -5
# 实例： 显示第3-5行和所有列数据
awk ' NR>=3 && NR <=5 {print $0}' TAIR10_GFF3_genes.gff
# 显示列数,-F指定分隔符，此处假定是table键分隔，默认空格键
awk -F "\t" '{print NF; exit}' some_data.bed
# 分隔符统一为一个空格，接下来将分隔符统一修改为逗号“，”
awk 'BEGIN{ FS=" ";OFS="," }{ print $1,$2,$3,$4 }' file2.txt > file3.txt
# link 
http://tiramisutes.github.io/2016/08/11/awk-forward.html
http://man.linuxde.net/awk
```

# Vim

```shell
https://draapho.github.io/2016/10/01/1604-CheatSheet-vim/
https://external-preview.redd.it/iigrixvxp5aYN9ox7Gr1dfI_rhLRotWlLsCafjJqjEQ.png?auto=webp&s=1594ddc17408cb9186a73c2a6d1a1bf1e00769dd
http://www.czhzero.com/2017/01/17/vim-key-shortcut/
# 
：help
# 启动
vim -c cmd file: 在打开文件前，先执行指定的命令；
vim -r file: 恢复上次异常退出的文件；
vim -R file: 以只读的方式打开文件，但可以强制保存；
vim -M file: 以只读的方式打开文件，不可以强制保存；
vim -y num file: 将编辑窗口的大小设为num行；
vim + file: 从文件的末尾开始；
vim +num file: 从第num行开始；
vim +/string file: 打开file，并将光标停留在第一个找到的string上。
vim --remote file: 用已有的vim进程打开指定的文件。 如果你不想启用多个vim会话，这个很有用。但要注意， 如果你用vim，会寻找名叫VIM的服务器；如果你已经有一个gvim在运行了， 你可以用gvim --remote file在已有的gvim中打开文件。
# 命令模式
ctrl+[
esc
vim test.c
# edit
i -- insert
ea -- end append
# quit with saving
:q               不保存退出
:q!              没权限情况下，强制不保存退出
:wq              保存修改并退出
:wq!             强行保存退出，(文件所属者科忽略文件的制度属性)
:w               保存不退出
:x or :exit      类似:wq，但是只有在有更改的情况下才保存
:qa              退出所有文件
:wqa             保存所有文件，退出
:w new_filename  另存为指定文件
ZZ               正常模式下，保存修改并退出(前面没有冒号)
ZQ               正常模式下，不保存退出
# 恢复
u                         取消上一步操作(最多连续取消500次)
Ctrl+r                    恢复上一步被撤销的操作
U                         撤销当前一行的操作
:e!                       返回上次保存后的状态
# 删除多行
# delete 1st to 100th lines
：set nu
1,100d
# delete
D        删除当前光标所在位置到某一行的结尾
d$       删除当前光标所在位置到某一行的结尾
dd       删除当前所在行
5dd     删除从当前行至其后的5行内容
dL       删除当前位置到屏幕上最后一行的内容
dH       删除当前位置到屏幕上第一行的内容
dG       删除当前位置到工作缓存区结尾的内容
d1G     删除当前位置到工作缓存区开始的内容
# delete commonts
1. Ctrl + V
x or d (如果是“//”注释，那需要执行两次该操作，如果是“#”注释，一次即可)
2. dd (删除光标所在的行)
# COPY (执行该操作首先需要进入视觉模式及ctrl+v, yank(提起))
yy复制游标所在行整行。或大写一个Y。
2yy或y2y复制两行.
y^复制至行首，或y0。不含游标所在处字元。
y$复制至行尾。含游标所在处字元。
yw复制一个word。
y2w复制两个字（单词）。
yG复制至档尾。
y1G复制至档首。
# CUT
Ctrl+v +dd
# PASTE
p or P
# search (n or N -> next or previous)
1. 光标停留在想要查找的单词的任意一个字母上面， 然后输入Shift + *
2. yw拷贝该单词， 然后输入 / (Ctrl + R) 0 （即 /”0）
3. / or ? 
# highlight search
Shift+*
set hlsearch
nohlsearch
# line number
set number
```

#  Emacs

```C++
http://ergoemacs.org/emacs/emacs_basics.html
https://www.gnu.org/
https://blog.csdn.net/redguardtoo/article/details/7222501
```

# Regular

```shell
# ^代表行首，$代表行尾。 ^$是空行的意思
# -v 表示反向选择
# 查询/etc/rsyslog.conf文件，但是不包含空行和注释行
grep -v '^$' /etc/rsyslog.conf | grep -v '^#'
```



# C

```c
gcc test.c
编译链接
gcc -o hello hello.c
```

# C++

```c++
// create file
vim test.cpp
touch test.cpp
// compile 
g++ test.cpp
// check version
g++ -v
// commonts
1. //
2. #
3. /*
   */
// set name for output file
g++ test.cpp -o output
// execute/run
./output
```

## For

```c++
// 
int main(){
    ...
    reture 0;
}
```



## wile

```C++
int count = 0;
while (!ireader.eof){
    reader >> item;
    if (item == "orange"){
        count ++;
    }
}
cout << count << " instancese of Orange found" << endl;
```

## read and write

```C++
#include <fstream>
//to read
1.
ifstream reader("file1.txt");
2.
ifstream reader;
reader.open("file.txt");
int a, b;
reader >> a >> b;
reader.close()
//to write
1.
ofstream writer("file1.txt");
2.
ofstream myfile;
myfile.open("newfile.txt")
writer << "some text" <<endl;
writer.close()
// if 
if(writer.is_open())
{
    ...
    cout << "Writing was success" << endl;
}
else
{
    cout << "Error" << endl;
}
2. 
if (reader.fail()){
    ...
}
```



# Server

### infor

```shell
#storage, memory,core 
cat /proc/cpuinfo |grep process |wc 
free -g
df -h
```



### upload/download

```shell
# To upload
#upload file
# SFTP
go to Terminal > Shell > New remote connection 
put /Users/liu/Desktop/Manuscript.docx /home/u2793/LiuGuojun/WGS_data
put /Volumes/LIU_WD/01Data/01PID/SHCHVD26042018_2.fastq.gz /home/u2793/LiuGuojun/WGS_data/PID_patient_1/copy
scp ~/local/file user@remote:~/file
scp -r 本地目录 linux用户名@ip地址:目标路径
man scp (查看文件)
ftp
```



### nohup

```shell
# quit
exit
# Quit the server (keep process running as background)
^a + ^d (control+z)
# pause the process
^z
# see the stoped file
jobs
# restart
fg %1 (%1 corresponds to the [1] with the jobs command)
fg %[number]
bg 
kill %[number]
kill -[signal] %[number]
disown %[number]
# running in background
nohup tar -czf iso.tar.gz Templates/* &
nohup ping google.com &
# process isn't killed when the terminal closes
disown  -h  %1
```

### Screen

```shell
sudo apt-get install screen
# go into screen
screen
screen -S bioinformatics (creating name)
screen -d -m python counter.py
# help
ctrl A > ?
# detach（退出）
ctrl A > Dscreen -S bioinfo
# show detach screen
screen -ls
# reatach 
screen -r # if a single screen
screen -d -r 10892.bioin # for attached screen
screen -r (number) # for detached screen
sudo screen -r
# kill screen
screen -X -S 12108 quit
ctrl a > k (kill screen)
ctrl a > \
# check if it is screen session
ps -e | grep 123123 (number)
# create window
ctrl a > c
# split screent
ctrl a + S (horizental window)
ctrl a + : + split (horizental window)
ctrl a + |
ctrl a + : + screen
ctrl a > : > split -v
ctrl a > : > focus
ctrl a > Q (quit screen)
ctrl a > : > only
ctrl a + H (create log file)
ctrl a + h (screenshot) 
# previous page , next page
ctrl + a + p
ctrl + a + n
ctrl + a + 1
ctrl + a + 2
ctrl + a + ""
# quit split window
ctrl a > X
ctrl a > : remove
# retrun
^A + N
# exit screen
exit
# title
ctrl a > A


```

### Tmux

```python
https://tmuxcheatsheet.com/
https://reishin.me/tmux/
https://zombie110year.top/2019/tmux-%E7%AE%80%E6%98%8E%E6%95%99%E7%A8%8B/
# basic
复制模式
Ctrl+b [ 空格标记复制开始，回车结束复制。
Ctrl+b ]   粘贴最后一个缓冲区内容
Ctrl+b =  选择性粘贴缓冲区
Ctrl+b :list-buffer  列出缓冲区目标
Ctrl+b :show-buffer 查看缓冲区内容
Ctrl+b :set mode-keys vi vi模式
Ctrl+b t 显示时间
# help
ctrl+b ?
tmux -V
Ctrl+b :list-commands
Ctrl+b :list-keys)
# To create a new named session
tmux new -s session_name
# check the session list
tmux ls
tmux list-sessions
Ctrl+b s (show all sessions)
# detach from the Tmux session
ctrl+b d
# Reattach
tmux a/at/attach/attach-session
tmux a/at/attach/attach-session -t mysession
# switch between sessions
Ctrl + b (
Ctrl + b )
# rename
tmux rename <newname>
ctrl+b $
# kill session
exit
tmux kill-ses -t mysession
tmux kill-session -t mysession
Ctrl+b :kill-server # kill all session
# kill sessions
tmux kill-session -a (kill all sessions but the current)
tmux kill-session -a -t mysession (kill all sessions but mysession)
# display name
tmux display-message -p '#S'
# split Panel
Ctrl+b c Create a new window (with shell)
Ctrl+b , Rename the current window
Ctrl+b % Split current pane horizontally into two panes
Ctrl+b " Split current pane vertically into two panes
Ctrl+b o Go to the next pane
Ctrl+b ; Toggle between the current and previous pane
Ctrl+b x Close the current pane
Ctrl+b z 最大化或者恢复面板
Ctrl+b <space>
Ctrl+b :resize-pane -U #向上调整
Ctrl+b :resize-pane -D #向下调整
  Ctrl+b :resize-pane -D 50
Ctrl+b :resize-pane -L #向左调整
Ctrl+b :resize-pane -R #向右调整
Ctrl+b ! 移动pane至window
Ctrl+b arrow key
# 移动pane合并至某个window
Ctrl+b :join-pane -t $window_name
Ctrl+b q # 显示pane编号
Ctrl+b Ctrl+o #按顺序移动pane位置
# window
# 切换到下一个窗口（n=next）
Ctrl+b c
Ctrl+b ,
# 切换窗口
Ctrl+b n
Ctrl+b p
Ctrl+b <num>  (切换到某个编号的窗口)
Ctrl+b 0 Switch to window 0 (by number )
Ctrl+b f # 如果窗口数量超过 9 个，通过窗口名称查找（f=find）
Ctrl+b w # Choose window from a list 如果窗口数量超过 9 个，通过窗口列表查找（w=window）
Ctrl+b l 在相邻的两个window里切换
# 关闭一个窗口
exit 或者 prefix &
Ctrl+b & Close current window
# 命令模式
Ctrl+b :
new-window -n console # 创建新窗口
new-window -n processes "top"`。 # 创建新窗口并执行命令

```



### Slurm

```R
# Useful links
http://parallel.uran.ru/book/export/html/313
https://slurm.schedmd.com/cpu_management.html
# check
squeue -u u2793
squeue --user=`whoami`
squeue --long
scontrol show job 174457
squeue --states=RUNNING
# clustering information
sinfo -s
srun --version
# detailed information about nodes, sections, tasks
scontrol show node tesla34
scontrol show partition
# run
sbatch bwa.sh
sbatch -c 12 -n 1 bwa.sh
sbatch -c 4 -n 3 bwa.sh
sbatch -N1 -n1 --mem-per-cpu=100m -t00:05:00 --qos=test hostname.sh
# cancel jobs
scancel 1897429 189723
man sbatch
# apply compute node
salloc -N 1 --exclusive -p work
ssh ...(nodelist/reason)
#
#SBATCH -C tesla  # choose tesla GPU
#SBATCH --gres=gpu:1 # use 1 GPU per node
#SBATCH --ntasks-per-node=12 #Set number of cores on the compute node.
#!/bin/sh
#SBATCH --partition=8CPUNodes # Name of the Slurm partition used
#SBATCH --job-name=STAR # Job Name
#SBATCH --nodes=5 # Select node
# 如果一个代码里跑四个程序。
#SBATCH --cpus-per-task=4 # 4 CPU allocation per Task
#SBATCH --time=20:00:00
# 内存设的越高越好，如果设的太高它不会被立即执行。
#SBATCH --mem=100GB
#SBATCH --wordir="/./././."
#SBATCH --error=<Error File Name>
#SBATCH --output=<Output File Name>

##########################
# example
#!/bin/bash
#SBATCH --job-name=Align
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --array=1-6
#SBATCH --time=20:00:00 
#SBATCH --output=bwa.%j.out
#SBATCH --mail-type=FAIL      # NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=gjliu0325@gmail.com
          
bwa mem -t 12 -M /home/u2793/_scratch/reference/index/bwa/hg19/bwa SHCHVD26042018_1.fastq.gz SHCHVD26042018_2.fastq.gz > bwa.sam
```

# Github

```shell
# 生成ssh key
ssh-keygen -t rsa -C "gjliu0325@gmail.com"
然后设置密码，我设置了，跟github账户差不多。最后少一个符号。
# help
which git
git --help
git status
# 查看远程主机的详细信息
git remote show origin
# create a new repository
echo "# Bash" >> README.md
git init
git add README.md
# 本地向远程github仓库提交文件(分三步)
1. 向本地stage增加文件,点号可以换成具体文件的名称（支持文件夹、通配符等）
git add . 
如果想撤销,使用git reset .（点号可以换成具体文件的名称（支持文件夹、通配符等））
或者使用git rm --cached <added_file_to_undo>
2. 向本地repos提交
git commit -m "提交日志"
3. 向远程github提交
git push -u origin master
如果本地某些文件不是最新的，可能需要先执行git pull更新一下（可能有冲突，需要自己手动合并一下，并填写合并日志）
# commit 
git commit
git commit -m "first commit"
# create a new repository on the command line
echo "# mynote" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin git@github.com:6guojun/mynote.git
git push -u origin master
# push an existing repository
git remote add origin https://github.com/guojun-liu/Bash.git
git remote add origin git@github.com:6guojun/mynote.git
git push -u origin master
# upload
git push -u origin master
# 命令用于添加远程主机
git remote add <主机名> <网址>
# 删除远程主机
git remote rm <主机名>
# 远程主机的改名
git remote rename <原主机名> <新主机名>
# 复制一个 Git 仓库
git clone https://github.com/guojun-liu/Bash.git
git clone <版本库的网址> <本地目录名>
git clone -o jQuery https://github.com/jquery/jquery.git # -0 设置远程主机名，默认是origin
# 远程主机的更新(commit)取回本地,合并远程分支
git fetch origin
git merge origin/master
#
git remote -v
#
ssh-keygen -t rsa -C
ssh -T git@github.comgit config --global user.email "jmzeng1314@163.com"
# 把所有文件加入到索引（不想把所有文件加入，可以用gitignore或add 具体文件)
git add - this stages your changes for committing
git add -A
git commit - this commits your staged changes locally
git push - this pushes your committed changes to a remotegit pullgit status
# 创建一个新的分支
git checkout -b newBrach origin/mastergit log
git pull 
# cancel
git reset
# delete
git rm
git rebase origin/master
# 取回远程主机某个分支的更新，再与本地的指定分支合并
git pull origin next:master
# 查看所有远程分支
git branch -a
# 命令指定master分支追踪origin/next分支
git branch --set-upstream master origin/next
git tag
# 准备删除本地和Git里的
git commit -a -m
//------------------------------常见错误-----------------------------------
1.$ git remote add origin git@github.com:WadeLeng/hello-world.git
错误提示：fatal: remote origin already exists.
解决办法：$ git remote rm origin
然后在执行：$ git remote add origin git@github.com:WadeLeng/hello-world.git 就不会报错误了
2. $ git push origin master
错误提示：error:failed to push som refs to
解决办法：$ git pull origin master //先把远程服务器github上面的文件拉先来，再push 上去。
# link
http://www.bio-info-trainee.com/2477.html
https://github.com/codepath/ios_guides/wiki/Using-Git-with-Terminalgit remote add origin https://github.com/yourUserName/yourRepoName.git
http://www.ruanyifeng.com/blog/2014/06/git_remote.html

```

