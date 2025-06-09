#install.packages("pheatmap")


library(pheatmap)         #???ð?
setwd("D:\\cell-death-iPF\\5数据集合并分析整套\\lasso建模验证ROC")      #???ù???Ŀ¼

#???????????ߵĺ???
bioRiskPlot=function(inputFile=null, project=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
	rt=rt[order(rt$riskScore),]      #???ղ??˵ķ??յ÷ֶ???Ʒ????????
		
	#???Ʒ???????
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)",
		 ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		
	#????????״̬ͼ
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)",
		 ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
	
	#??????ͼע?͵???ɫ
	ann_colors=list()
	bioCol=c("blue", "red")
	names(bioCol)=c("low", "high")
	ann_colors[["Risk"]]=bioCol

	#???Ʒ?????ͼ
	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(Risk=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation,
		     annotation_colors = ann_colors, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("blue",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}

#???ú???,????tarin??????????
bioRiskPlot(inputFile="70866trainRisk.txt", project="70866train")
#????test??????????
bioRiskPlot(inputFile="70867risk-real.txt", project="70867test")
#??????????Ʒ????????
bioRiskPlot(inputFile="70866TestRisk-real.txt", project="70866test")


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

