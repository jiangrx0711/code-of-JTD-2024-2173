######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("pheatmap")

library(pheatmap)
setwd("/Users/limeng/Desktop/IPF-ICD-Final/模型构建28/模型30最终/风险热图")                   #设置工作目录

#绘制train组风险热图
rt=read.table("riskTrain.txt",sep="\t",header=T,row.names=1,check.names=F)      #读取train输入文件
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
rt1=log2(rt1+0.01)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmapTrain.pdf",width = 12,height = 5)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("deepskyblue2", "cadetblue", "azure","darkgoldenrod2","brown3"))(50) )
dev.off()

#绘制test组风险热图
rt=read.table("riskTest.txt",sep="\t",header=T,row.names=1,check.names=F)      #读取test输入文件
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
rt1=log2(rt1+0.01)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmapTest.pdf",width = 12,height = 5)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("deepskyblue2", "cadetblue", "azure","darkgoldenrod2","brown3"))(50) )
dev.off()

######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056