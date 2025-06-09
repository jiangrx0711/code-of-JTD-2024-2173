
library(venn)                  
outFile="intersectGenes.txt"        
setwd("/Users/limeng/Desktop/IPF-ICD/Venn")    #
geneList=list()

rt=read.table("diffGeneExp.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,1])              
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)             
geneList[["DEGs"]]=uniqGene

rt=read.table("tcga.uniCox.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,1])               
geneNames=gsub("^ | $","",geneNames)    
uniqGene=unique(geneNames)             
geneList[["Prognostic genes"]]=uniqGene


mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
