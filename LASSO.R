library(tidyverse)
setwd("D:\\BaiduNetdiskDownload\\LUAD-硫化\\lasso")
geneCoef <- data.table::fread('geneCoef.txt')
exp <- data.table::fread('30219+31210-time-exp.txt')

RS=c()
for (i in 1:309) {
  x <- as.numeric(geneCoef[1,2]*exp[i,4])+as.numeric(geneCoef[2,2]*exp[i,5])+
    as.numeric(geneCoef[3,2]*exp[i,6])+as.numeric(geneCoef[4,2]*exp[i,7])
  RS[i] <- x
}

exp$RS <- RS
exp$group <- ifelse(exp$RS>median(exp$RS),'High','Low')
table(exp$group)
write.csv(exp,file = '30219+31210-time-exp.txt')
