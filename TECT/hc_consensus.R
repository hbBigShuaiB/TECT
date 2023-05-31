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

# #————————————————————————
# inc<-read.csv('/Users/yzf/PycharmProjects/TECT/test/mycount.csv',header=TRUE,row.names=1)
# inc<-t(inc)  # 转为基因x样本
# inc<-as.data.frame(inc)
# sum_exon<-inc['sum_exon',]
# sum_map<-inc['sum_map',]
# inc<-inc[-grep('sum_exon',rownames(inc)),]
# inc<-inc[-grep('sum_map',rownames(inc)),]
# #————————————————————————
inc<-read.csv(myArgs[1],header=TRUE,row.names=1)
sum_exon<-inc[,'sum_exon']
sum_map<-inc[,'sum_map']
inc<-inc[,-grep('sum_exon',colnames(inc))]
inc<-inc[,-grep('sum_map',colnames(inc))]
inc<-t(inc)  # 转为基因x样本

# sum_exon<-inc['sum_exon',]
# sum_map<-inc['sum_map',]
# inc<-inc[-grep('sum_exon',rownames(inc)),]
# inc<-inc[-grep('sum_map',rownames(inc)),]

x<-DGEList(inc) #给入DGE
gsm_sample <- substring(colnames(x), 1, nchar(colnames(x))-16) # 得到gsm+样本名
# 给上样本分组和缩减细胞gsm
group<-substring(gsm_sample,12,nchar(gsm_sample))
x$samples$group<-group
samplenames <- substring(colnames(x), 1, 10)
colnames(x)<-samplenames

#给上celltype
celltype<-read.table(myArgs[2],header=TRUE,row.names = 1)
# celltype<-read.table('../test/cell_type.txt',header=TRUE,row.names = 1)
x$samples$celltype<-c('null')
for(i in rownames(x$samples)){
  x$samples[i,]$celltype<-celltype[i,]
}
# 给上condition
condition<-read.table(myArgs[3],header=TRUE,row.names=1)
# condition<-read.table('../test/condition.txt',header=TRUE,row.names=1)
x$samples$condition<-c('null')
for(i in rownames(x$samples)){
  x$samples[i,]$condition<-condition[i,]
}
condition<-x$samples$condition
# 给上sum_exon,sum_map
x$samples$sum_exon<-sum_exon
x$samples$sum_map<-sum_map

colnames(inc)<-substring(colnames(inc),1,10)





print('开始过滤')
# 过滤
# 某细胞表达的基因种类数小于600就过滤,去掉了73个细胞 mgt
for(i in colnames(x$counts)){
  if(length(which((x$counts[,i]!=0)))<mgt){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

print('第二步')
#过滤掉总count数过少或过多的细胞 minc maxc
#x$counts<-x$counts[,colSums(x$counts)<1000000] 
#x$counts<-x$counts[,colSums(x$counts)>100000] 
## 太少的过滤掉15个
for(i in colnames(x$counts)){
  if(sum(x$counts[,i])<minc){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

print('第三步')
## 太多的过滤掉18个
for(i in colnames(x$counts)){
  if(sum(x$counts[,i])>maxc){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

print('第四步')
#write.csv(x$counts,'tes.csv') mgc
#过滤掉表达量低的基因
b<-c()
for(i in 1:nrow(x$counts)){
  if(sum(x$counts[i,])<mgc){
    b<-append(b,i)
  }
}
if(!is.null(b)){
  x$counts<-x$counts[-b,]
}
print('过滤结束')



cell_types<-unique(x$samples$celltype)
consensus_cluster <- function(cell_type) {
  this_title<-paste0("result/hc/",cell_type)
  print(this_title)
  some_type<-inc[,which(x$samples$celltype==cell_type)]
  mads=apply(some_type,1,mad)
  some_type<-some_type[rev(order(mads))[1:500],]
  some_type<-sweep(some_type,1,apply(some_type,1,median))
  ConsensusClusterPlus(some_type,maxK=6,reps=50,pItem=0.8,pFeature=1,title=cell_type,clusterAlg='hc',distance='pearson',plot='png')
}

for (cell_type in cell_types){
  consensus_cluster(cell_type)
}





# # 整体用重复序列进行一次聚类
# 
# d<-inc
# ## 计算绝对中位数（Median Absolute Deviation）
# mads=apply(d,1,mad)   # 每行——每个基因一个值
# d=d[rev(order(mads))[1:500],]   # order对mads进行从小到大排序，rev进行逆转——从大到小，然后取前500个大的数值
# d = sweep(d,1, apply(d,1,median,na.rm=T)) #  每一行的每个数减去这一行的中位数
# ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature = 1,title='210714.hc',clusterAlg='hc',distance='pearson',plot='png',writeTable=T)











# # 用beta细胞进行聚类
# some_type<-inc[,which(x$samples$celltype==cell_type)]
# mads=apply(some_type,1,mad)
# some_type<-some_type[rev(order(mads))[1:500],]
# some_type<-sweep(some_type,1,apply(some_type,1,median))
# ConsensusClusterPlus(some_type,maxK=6,reps=50,pItem=0.8,pFeature=1,title=cell_type,slusterAlg='hc',distance='pearson',plot='png',writeTable=T)



















