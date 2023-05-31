myArgs<-commandArgs(trailingOnly = TRUE)

mgt<-as.integer(myArgs[4])
minc<-as.integer(myArgs[5])
maxc<-as.integer(myArgs[6])
mgc<-as.integer(myArgs[7])

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ConsensusClusterPlus","limma","Glimma","edgeR","ggpubr","factoextra"))

library(limma)
library(Glimma)
library(edgeR)
library(ggpubr)
library(factoextra)
library(ConsensusClusterPlus)

# options(ggrepel.max.overlaps = Inf)
options(digits=3)

inc<-read.table(myArgs[1],header=TRUE,row.names=1)  # 1346*1600
inc<-t(inc)  # 转为基因x样本

sum_exon<-inc['sum_exon',]
sum_map<-inc['sum_map',]
inc<-inc[-grep('sum_exon',rownames(inc)),]
inc<-inc[-grep('sum_map',rownames(inc)),]

x<-DGEList(inc) #给入DGE
gsm_sample <- substring(colnames(x), 1, nchar(colnames(x))-16) # 得到gsm+样本名
# 给上样本分组和缩减细胞gsm
group<-substring(gsm_sample,12,nchar(gsm_sample))
x$samples$group<-group
samplenames <- substring(colnames(x), 1, 10)
colnames(x)<-samplenames

#给上celltype
celltype<-read.table(myArgs[2],header=TRUE,row.names = 1)
x$samples$celltype<-c('null')
for(i in rownames(x$samples)){
  x$samples[i,]$celltype<-celltype[i,]
}
# 给上condition
condition<-read.table(myArgs[3],header=TRUE,row.names=1)
x$samples$condition<-c('null')
for(i in rownames(x$samples)){
  x$samples[i,]$condition<-condition[i,]
}
condition<-x$samples$condition
# 给上sum_exon,sum_map
x$samples$sum_exon<-sum_exon
x$samples$sum_map<-sum_map
colnames(inc)<-substring(colnames(inc),1,10)



# 过滤
# 某细胞表达的基因种类数小于600就过滤,去掉了73个细胞
for(i in colnames(x$counts)){
  if(length(which((x$counts[,i]!=0)))<mgt){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

#过滤掉总count数过少或过多的细胞
#x$counts<-x$counts[,colSums(x$counts)<1000000]
#x$counts<-x$counts[,colSums(x$counts)>100000]
## 太少的过滤掉15个
for(i in colnames(x$counts)){
  if(sum(x$counts[,i])<minc){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

## 太多的过滤掉18个
for(i in colnames(x$counts)){
  if(sum(x$counts[,i])>maxc){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

#write.csv(x$counts,'tes.csv')
#过滤掉表达量低的基因
b<-c()
for(i in 1:nrow(x$counts)){
  if(sum(x$counts[i,])<mgc){
    b<-append(b,i)
  }
}
x$counts<-x$counts[-b,]

cell_types<-unique(x$samples$celltype)
consensus_cluster <- function(cell_type) {
  this_title<-paste0("result/hc/",cell_type,"/")
  print(this_title)
  some_type<-inc[,which(x$samples$celltype==cell_type)]
  mads=apply(some_type,1,mad)
  some_type<-some_type[rev(order(mads))[1:500],]
  some_type<-sweep(some_type,1,apply(some_type,1,median))
  ConsensusClusterPlus(some_type,maxK=6,reps=50,pItem=0.8,pFeature=1,title=this_title,clusterAlg='km',distance='pearson',plot='png')
}

for (cell_type in cell_types){
  consensus_cluster(cell_type)
}
