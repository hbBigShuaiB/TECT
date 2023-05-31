library(limma)
library(Glimma)
library(edgeR)
library(ggpubr)
library(factoextra)

# TECT/kmeans.R test/mycount.csv test/cell_type.txt test/condition.txt 600 10000 100000 10000 2      

myArgs<-commandArgs(trailingOnly = TRUE)
mgt<-as.numeric(myArgs[4])
minc<-as.numeric(myArgs[5])
maxc<-as.numeric(myArgs[6])
mgc<-as.numeric(myArgs[7])
options(digits=3)

cm<-read.csv(myArgs[1],header=TRUE,row.names=1)
cm<-t(cm)

sum_exon<-cm['sum_exon',]
sum_map<-cm['sum_map',]
cm<-cm[-grep('sum_exon',rownames(cm)),]
cm<-cm[-grep('sum_map',rownames(cm)),]

x<-DGEList(cm) #给入DGE
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

# 过滤
print('开始过滤')
# 某细胞表达的基因种类数小于600就过滤,去掉了73个细胞
for(i in colnames(x$counts)){
  if(length(which((x$counts[,i]!=0)))<600){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

#过滤掉总count数过少或过多的细胞
#x$counts<-x$counts[,colSums(x$counts)<1000000]
#x$counts<-x$counts[,colSums(x$counts)>100000]

for(i in colnames(x$counts)){
  if(sum(x$counts[,i])<100000){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}
for(i in colnames(x$counts)){
  if(sum(x$counts[,i])>1000000){
    x$counts<-x$counts[,-grep(i,colnames(x$counts))]
    x$samples<-x$samples[-grep(i,rownames(x$samples)),]
  }
}

#write.csv(x$counts,'tes.csv')
#过滤掉表达量低的基因
b<-c()
for(i in 1:nrow(x$counts)){
  if(sum(x$counts[i,])<10000){
    b<-append(b,i)
  }
}
x$counts<-x$counts[-b,]

print('过滤结束')
cpm<-cpm(x)
cpm<-t(cpm)

df <- scale(cpm)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 2, linetype = 2)
fviz_nbclust(df, kmeans, method ="silhouette")  # 轮廓系数

km_result <- kmeans(df, 2, nstart = 24)
print(km_result)
addcluster <- cbind(df, cluster = km_result$cluster)
table(addcluster[,'cluster'])
picture<-fviz_cluster(km_result, data = df,
                      geom='point',
                      palette = c("#2E9FDF", "#00AFBB"),
                      ellipse.type = "euclid",
                      #star.plot = TRUE, 
                      repel = TRUE,
                      ggtheme = theme_minimal()
)

ggsave(picture,filename='kmeans.pdf',width=12,height=9)

