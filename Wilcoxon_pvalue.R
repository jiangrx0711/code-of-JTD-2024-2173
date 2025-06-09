#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

#引用包
library(limma)
library(pheatmap)

inputFile="normalize.txt"        #输入文件
pFilter=0.05                     #pvalue临界值
logFCfilter=0                    #logFC临界值
conFile="sample1.txt"            #对照组样品
treatFile="sample2.txt"          #实验组样品
setwd("/Users/limeng/Desktop/IPF-ICD-Final/2.差异分析")     #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

#读取样品信息
sample1=read.table(conFile,sep="\t",header=F,check.names=F)
sample2=read.table(treatFile,sep="\t",header=F,check.names=F)
conData=data[,as.vector(sample1[,1])]
treatData=data[,as.vector(sample2[,1])]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)
type=c(rep(1,conNum),rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
	if(sd(data[i,1:conNum])==0){
		data[i,1]=0.001
	}
	if(sd(data[i,(conNum+1):(conNum+treatNum)])==0){
		data[i,(conNum+1)]=0.001
	}
	geneName=i
	rt=data.frame(expression=data[i,],type=type)
	wilcoxTest<-wilcox.test(expression ~ type, data=rt)
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	pvalue=wilcoxTest$p.value
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,logFC=logFC,conMean=conGeneMeans,treatMean=treatGeneMeans,pValue=pvalue))
	}
}

#输出所有基因的差异情况
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$pValue))<pFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
write.table(outDiff,file="diff.txt",sep="\t",row.names=F,quote=F)

#差异基因的表达文件
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

#绘制差异基因热图
geneNum=50
diffSig=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(diffSig[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2) ){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=7,width=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
