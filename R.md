#1. Basic

```R
查看文件**
查看属性：str()
Language in terminal:
      search for "Renviron" in R folder, and add "LANGUAGE=en"
查看属性的数据：head(...$edgeData)
查看文件 dir()
查看单列数据：summary(as.factor(sex))
     dim()
     help.search(“t test”)
    ?read.csv
     view()
     mi_expr[1:100,1:100] %>% View()
     head()
     help(mean)
     is.numeric(rt)
     class(rt)
     fix(rt)
# Error: vector memory exhausted (limit reached?)
rm(list=ls())
ls()
rm(lissNPC,lissNPC.sta,lissNPC1.ady,lissNPC1.ref)
#
mathmatics
     demo(plotmath)          
tmp <- a[1:20]
**工作环境**
setwd(“”)
中断工作再继续工作：save(fpkm,file=‘fpkm_56breastcancer.Rdata ’)
load(‘fpkm_56breastcancer.Rdata’)
**画图**
cex: character expansion
bty: box type 
#读写数据
   rt=read.csv(“…”,header=T, sep=“,”,row.names=1, check.names=F)
   header是第一行变列名，
   rt=read.table(“…”,header=T, sep=“\t”, row.names=1, check.names=F)
  写文件：write.csv(a,file = “gene_expr.csv”,na=“”,row.name=FALSE, sep=“,”)
   rt  <- get(load(paste0 (“KIRC_CNV_results.rda")))
   save.image("BLCA.Rdata")
   load("BLCA.Rdata")
library(“readxl")
data <- read_excel("Normaltissue.xlsx", sheet=1, col_names=TRUE, col_types= NULL,na = "", skip = 0)
my_data <- read_excel("my_file.xlsx")

# if "package 'xxx' is not available then use the following method
source("http://bioconductor.org/biocLite.R")
biocLite("")
#Uninstall
remove.packages("")

```

### shortcut

```R
Shift + command + M  "%>%"
alt + - "<-"
alt + shift + K "manu"
```

### Upgrade

```R
# Check version 
R.Version()
# Update
install.packages('devtools') #assuming it is not already installed
library(devtools)
install_github('andreacirilloac/updateR')
library(updateR)
updateR(admin_password = 'lgj19890707')
#  transport R packages to the upgrade R
/Library/Frameworks/R.framework/Versions/x.xx/Resources/library
update.packages(checkBuilt=TRUE)
version
packageStatus()
# update all packages
update.packages(repos = "https://mirrors.ustc.edu.cn/CRAN/",ask='graphics',checkBuilt=TRUE)
# 删除旧版本的R包
/Library/Frameworks/R.framework/
```



# 2. Data frame

### NA value and delete 0, sapce

```R
# 用0替换NA值
rt[is.na(rt)] <- 0
# 数据插补：用中位数补全NA
result=data.imputation(Data,fun="median")
#
exprSet2[complete.cases(exprSet2),]
sum(complete.cases(exprSet2))
sum(!complete.cases(exprSet2))
# function
delete.na <- function(DF, n = 0){
  DF[rowSums(is.na(DF)) <= n,]
}
delete.na(exprSet2,1)
# 
library(tidyr)
exprSet2 %>% drop_na()
exprSet2 %>% drop_na(5:50)
#
row.has.na <- apply(final, 1, function(x){any(is.na(x))})
sum(row.has.na)
final.filtered <- final[!row.has.na,]
#
dplyr::filter(df,!is.na(columnname))
# delete empty rows, not NAs, you can do:
data[!apply(data == "", 1, all),]
#To remove both (NAs and empty):
data <- data[!apply(is.na(data) | data == "", 1, all),]
# delete 0
# keep rows have 0 less then 5 columns 
a <- mi_expr[apply(mi_expr[,-1] == 0, 1, sum) <= 5, ]
a <- mi_expr[apply(mi_expr[,-1],1,function(X) length(X[X>=50])>0 ),]
b <- df[apply(df[,-1], 1, function(x) !all(x==0)),]
# delete space
df$gene_id <- gsub('\\s+','',df$gene_id)
# select rows and columns
mRNA_lncRNA <- rt[,c(2:4,8:10)]
rt1 <- rt[c(2, 5:134),]
```

###Random data

```R
# column
df$status <– sample(0:1, 35, TRUE, prob = c(.7, .3))

# 建立矩阵
a <- matrix(1:20, nrow=5, ncol=4, byrow=T)
转换向量为矩阵
a <- 1:20
dim(a) <- c(5,4)
# Data frame 
a <- data.frame(replicate(5,sample(0:5,10,rep=TRUE)))
a <- a[order(a$X5),]
# order by column name
mRNA3 <- mRNA2[,order(colnames(mRNA2))]
#
school<-c("NYU", "BYU", "USC","AAA","BBB")
state<-c("NY","UT","CA","DF","AS")
measure<-c("MSAT","MSAT","GPA","ASD","ADA")
score<-c(1, 2, 2,4,4)
df<-data.frame(school,state, measure,score,  stringsAsFactors=FALSE)
#
df <- data.frame(id = rep(1:4, rep(2,4)),
                 visit = I(rep(c("Before","After"), 4)),
                 x = rnorm(4), y = runif(4))
#
dplyr::data_frame(
    ID = wakefield::id(n=10),
    Smokes = smokes(n=10),
    Sick = ifelse(Smokes, sample(5:10, 10, TRUE), sample(0:4, 10, TRUE)),
    Death = ifelse(Smokes, sample(0:1, 10, TRUE, prob = c(.2, .8)), sample(0:1, 10, TRUE, prob = c(.7, .3)))
)
# 
wdata = data.frame(
  sex = factor(rep(c("F", "M"), each=200)),
  weight = c(rnorm(200, 55), rnorm(200, 58))
)
head(wdata, 4)
#
df <- read.table(textConnection("
A1  A1  0.90
A1  B1  0.85
A1  C1  0.45
A1  D1  0.96
B1  B1  0.90
B1  C1  0.85
B1  D1  0.56
C1  C1  0.55
C1  D1  0.45
D1  D1  0.90"))
```

###String process

```R

substr
split: 将a文件按sample列的种类为标准进行分类，
      b = split(a, a$sample)  
      b[[1]]或b[[2]]
strsplit
trtrssplit

%>%: 管道函数，左件的值送给右件的表达式，并作为函数的第一个函数。
修改变量名：
 1.fix() 
 2.rename(dateframe,c(oldname=“newname”,oldname=“newname”,…))   
批量改变变量名：
       unlist: tranform list to vector
        lappy:
        修改列名：GSM1384316_ZR78B.txt
ZR78B
           names(fpkm) <- unlist(lapply(names(fpkm),function(x){
                            tmp <- strsplit(x,’_')[[1]][2]   #[2]留右边的
                            tmp <- strsplit(tmp,’\\.')[[1]][1]        }))
```

### Rowname & Columns names

```R
#Rownames
# OR4F5 | ENSG00000186092 | protein_coding to OR4F5
library(tidyr)
library(dplyr)
tmp <- tmp %>% 
  tidyr::separate(gene_id,into =c("gene_id","drop"),sep="\\|") %>%
  dplyr::select(-drop)
# TCGA.CU.A0YN.11A.11R.A10U.07 to TCGA.CU.A0YN
# gene_id 是第一列，不能是行名
rt$id <- unlist(lapply(rt$id,function(x){
  tmp <- substr(x,start = 1,stop = 12)
}))
sur$sample <- stringr::str_sub(sur$sample, start = 1 ,end= 12)
# TCGA.CU.A0YR to TCGA-CU-A0YR
rt$id <– gsub("\\.", "-", rt$id)
#Columns names
colnames(genomicmatrix) <- gsub("[.]", "_",colnames(genomicmatrix))
# set rownames
rownames(A)=A[,1]
# suppose first row is characters
colnames（a）=a[1,]
colnames(tmp3) <- as.character(unlist(tmp3[1,]))
# If not characters
tmp3[] <- lapply(tmp3, as.character)
# Row name to column
library(tibble)
library(dplyr)
tmp <- as.data.frame(tmp)
tmp <- tmp %>% rownames_to_column("ID")
head(tmp)
# Column name
# XGSE12333 to GSE12333
names(df) <- substring(names(df), 2)
# first column name
names(rt)[1] <- "V1"
# first row as column name
colnames(rt1) <- as.character(unlist(rt1[1,]))
rt1 = rt1[-1, ]
```

### List

```R
# combine two list
new <– c(a,b)
```



###list to dataframe

```r
#1
c2<- melt(c1,
          value.name = "entropy",
        # id = 1,
        # measure = c("","") 
          na.rm=T,
          variable.names("entropy")
          )
#2
df <- data.frame(matrix(unlist(c1), nrow=k, byrow=T),stringsAsFactors=TRUE)
df2 = do.call(rbind.data.frame, c1)
# to numeric
rownames(a) = a[,1]
a <- data.matrix(a)
a <- t(a)

```

###Add

```R
# add the element of two columns
df$x <- paste(df$n,df$s)
# add two column 
A <- data.frame(tra[,c(1,2)])
B <- data.frame(tra[,c(3,4)])
C<- rbind(A,B)
取列名/行名：
colnames(rt) = paste("Sample-name", 1:10, sep = "")#生成gene_1这样文   件
rownames(rt) = paste("Gene-name", 1:20, sep = "")
names(datExpr0) = blcaData$substanceBXH
names(datExpr0)=blcaData$trans
library(reshape2)        
fpkm <- dcast(a,formula = V2~V1)
melt:
计算单行键放新变量里：
\1. mydata$meanx <- (mydata$x1 + mydata$x2) /2
         \2. transform(mydata, sumx= x1 +x2,…,….)
更改列的类型：
Affairs$ynaffair[Affairs$affairs >0]<-1
        Affairs$ynaffair[Affairs$affairs ==0]<-0
数据框运算：
apply：returns vectors or list
mads=apply(d,1,mad)  #计算中位数差,  1 是行，2是列
sapply: returns only vector as output
lapply: returns only list as output
mapply: multivariate version of sapply
tapply: apply a function to **subsets** of a vector
rapply: recursive apply
vapply: similar to sapply, but has a pre-specified type of return value
选取某列：df <- subset(df, select = c(a,c))
删除某列：df <- subset(df, select = -c(a,c) )
标准化sweep
 d=sweep(d,1, apply(d,1,median,na.rm=T)) #按行对所有数据都减去各列的平均数
按某列排序后删除异常值：
LncRNA <- LncRNA[order(LncRNA[,1]),] 
转换数值型数据：
class(rt)
rt[] <- lapply(rt, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(rt, class)
# unite two columns
a$new <- paste(a$GSE51406,a$CVID,sep = "_")
```

### Delete

```R
# Delete columns by names
tmp3<- tmp2[,-which(names(tmp2) %in% c("value"))]
# Delete columns by number
去掉某列：rt<-rt[,-1] 
数据框中提取某些列：
# 留：1, 2 和5 列
keep <- seq(2, ncol(rt),3)
count_data <- cbind(rt[,1],rt[,keep])
LncRNA <- LncRNA[-(1:9),]
# delete normal sample
index_mRNA <- as.numeric(substr(colnames(mRNA),start = 14, stop = 15))
mRNA1 <-  mRNA[,which(index_mRNA  == 1)] 
```



###Match

```R
两矩阵提取相同部分：
c <- a[match(b$id, a$id),] # b为标准, id 是第一列名
phenoData[match(rownames(exprs), rownames(phenoData)),]#its first argument **in** its second.
# match columns and keep all duplicated elements
c = r2[r2$V1 %in% r1$V1,]  # r1为标准
phenoData = phenoData[rownames(phenoData) %**in**% rownames(expr),]# match for its left operand
合并矩阵：merge_data = merge(expr, phenoData, by=0, all.x=T) #0以行
对某行进行运算：
cast: 把重复的多行少列矩阵变成V1是行,V2是列的矩阵,
c <- intersect(a$SYMBOL,b$logFC)
# match multiple columns
new <- Reduce(intersect, list(a,b,c,d))
# match by column
mRNA2 <- mRNA1[,which(colnames(mRNA1) %in%  sur$sample)]

```



###Deduplicate elements:

```R
#Deduplicate
index <- duplicated(lncRNA$gene_id)
LncRNA <- lncRNA [!index,]
# Or
unique(a)
# remove unique and keep duplicates
data <- data.frame(id1 = c(1, 1, 1, 2, 2),
                   id2 = c(2, 2, 3, 2, 2))
data[data$id2 %in% data$id2[duplicated(data$id2)],]
#keep duplicates
library(data.table)
setkey(setDT(data),id2) # or with data.table 1.9.5+: setDT(dx,key="ID")
data[duplicated(data) |duplicated(data,fromLast=T)]
#keep duplicates
data[ ave(1:nrow(data), data$id2, FUN=length) > 5 , ]
#keep duplicates
data[ unlist(tapply(1:nrow(data), data$id2, function(x) if (length(x) > 1) x)) , ]
#keep duplicates
library(dplyr)
df %>% group_by(ID) %>% filter( n() > 1 )

```

### Transform (long & wide)

```R
# long to wide
m <- dcast(m,formula = Core_Gene~Target, value.var = "Distance")
```

### Converting data.frame 

```R
# to numerical matrix

```



# 3. Plot

```R
pdf("A.pdf")  dev.off()
text(5, 10.2, "..."/)
text(4, 8.4, "expression(hat(beta) == (X^t * X)^{-1} * X^t * y)", cex = .75)
mtext("«Latin-1 accented chars»: é è ø å æ", side = 3)
plot(1:10, 1:10, main = “text….”, sub = “….”)

```

#                          4. Loop

```R
# while 
i <- 1
sum <- 0
while(i <= 100){sum = sum + 1; i = i + 1}
print(sum)
# ifelse
mydata$x4 = ifelse(mydata$x2>150,1,0)
ifelse(mydata$x1<10 & mydata$x2>150,1,0) # "And" condition
ifelse(mydata$x1<10 | mydata$x2>150,1,0) # "Or" Condition
mydata$y = ifelse(mydata$x3 %in% c(“A”,”D”) ,mydata$x1*2,mydata$x1*3)
ifelse(is.na(x),1,0)
# Nested ifelse Statement
mydata$y = ifelse(mydata$x3 %in% c(“A”,”B”) ,mydata$x1*2,
                  ifelse(mydata$x3 %in% c(“C”,”D”), mydata$x1*3,
                         mydata$x1*4))
### for in vector or list 
for(i in seq(from = 1, to = 100, by = 1)) sum = sum + i 
print(sum)
# if
k = 100
if (k > 100){
  print(" > 100 ")
} else if (k < 100){
  print(" < 100 ")
} else {
  print(" =100 ")
}
# iterators
iter_one <- iter(1:10, checkFunc = function (x) x%%2==0, recycle = F)
nextElem(iter_one)
# foreach, traversing the vector matrix dataframe or iterator
library(foreach)
i_square <- foreach(i = 1:5) %do% i^2
i_square
# dplyr
library(dplyr)
x=c(1,NA,2,3)
if_else(x%%2==0, “Multiple of 2”, “Not a multiple of 2”, “Missing”)
# SQL grammar
df=data.frame(k=c(2,NA,3,4,5))
library(sqldf)
sqldf(
  "SELECT *,
  CASE WHEN (k%2)=0  THEN 'Multiple of 2'
  WHEN  k is NULL  THEN 'Missing'
  ELSE 'Not a multiple of 2'
  END AS T
  FROM df"
)

```

### Progress bar

```R
library(progress)
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = 100, clear = FALSE, width= 60)
for (i in 1:100) {
  pb$tick()
  Sys.sleep(1 / 100)
}
```

| matlab | python | R    |
| ------ | ------ | ---- |
|        |        |      |

### Multithreading

```R
Check core:

      detectCores()

      detectCores(logical = FALSE)

parallel: 主要针对lappy 函数

    stopCluster(cl)    

foreach: 针对for 函数

library(foreach)

library(doMC)

registerDoMC(cores=4)

a <- foreach(b = vec, .combine=rbind) %dopar% {

  data.frame(b, c = b^2)

}

for:

    GPU:

    gpuR
# library(doParallel)
# library(foreach)
# cores=detectCores()
# cl <- makeCluster(cores[1]-4) #not to overload your computer
# registerDoParallel(cl)
# foreach(s = 1:ncol(tmpE.m)) %do% {
#   pb$tick()
#   sr.lo[[s]] <- CompSR(tmpE.m[,s],pin.o$a)
#   print(s)
# }
# stopCluster(cl)
# stopImplicitCluster()

```

# R package

```R
library('devtools') # 开发 R 包黑魔法工具
create('~/somebm') # 建立 R 包的目录， somebm 就是你想要的包的名称
setwd('~/somebm') # 把工作目录放到 R 包中的目录，开发 R 包过程中始终推荐这样做。
dir() # 列出当前工作目录的文件和文件夹
```



