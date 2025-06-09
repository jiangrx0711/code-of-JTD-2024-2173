library(survival)
library(survminer)
library(timeROC)




setwd("D:\\乳酸LUAD\\roc") 



bioROC=function(inputFile=null,rocFile=null){
	
	rt=read.table(inputFile,header=T,sep="\t")

	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,2,3),ROC=TRUE)
	pdf(file=rocFile,width=5,height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

#bioROC(inputFile="tcgaRisk.txt",rocFile="tcga.ROC.pdf")
bioROC(inputFile="622-5risk (2).txt",rocFile="622-5risk (2)ROC.pdf")


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
