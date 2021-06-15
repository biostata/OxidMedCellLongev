

#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("survivalROC")

#引用包
library(survival)
library(caret)
library("glmnet")
library(survminer)
library("survivalROC")
# library(dplyr)

            #设置工作目录
setwd("C:\\Users\\lenovo\\Desktop")                 #设置工作目录

rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取输入文件
rt[,"futime"]=rt[,"futime"]/365                                               #生存时间单位改为年
# rt=rt %>% mutate_if(is.character, as.factor)
rt<- as.data.frame(rt, stringsAsFactors = F , na.rm = T)

#对分组进行循环，找出train和test都显著的分组
for(i in 1:1000){
	  #############对数据进行分组#############
		inTrain<-createDataPartition(y=rt[,2],p=0.5,list=F)
		train<-rt[inTrain,]
		test<-rt[-inTrain,]
		trainOut=cbind(id=row.names(train),train)
		testOut=cbind(id=row.names(test),test)
		
		#############单因素COX分析#############
		outTab=data.frame()
		pFilter=0.05
		sigGenes=c("futime","fustat")
		for(i in colnames(train[,3:ncol(train)])){
					 cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
					 coxSummary = summary(cox)
					 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
					 outTab=rbind(outTab,
					              cbind(id=i,
					              HR=coxSummary$conf.int[,"exp(coef)"],
					              HR.95L=coxSummary$conf.int[,"lower .95"],
					              HR.95H=coxSummary$conf.int[,"upper .95"],
					              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
					              )
					  if(coxP<pFilter){
					      sigGenes=c(sigGenes,i)
					  }
		}
		train=train[,sigGenes]
		test=test[,sigGenes]
	  uniSigExp=train[,sigGenes]
	  uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
	
	  #############lasso回归#############
		trainLasso=train
		trainLasso$futime[trainLasso$futime<=0]=0.003
		x=as.matrix(trainLasso[,c(3:ncol(trainLasso))])
		y=data.matrix(Surv(trainLasso$futime,trainLasso$fustat))
		fit <- glmnet(x, y, family = "cox", maxit = 1000)
		cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
		coef <- coef(fit, s = cvfit$lambda.min)
		index <- which(coef != 0)
		actCoef <- coef[index]
		lassoGene=row.names(coef)[index]
		lassoGene=c("futime","fustat",lassoGene)
		if(length(lassoGene)==2){
		   next
		}	
		train=train[,lassoGene]
		test=test[,lassoGene]
		lassoSigExp=train
	  lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
	
	  #############构建COX模型#############
	  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
	  multiCox=step(multiCox,direction = "both")
	  multiCoxSum=summary(multiCox)
		
		#输出模型相关信息
		outMultiTab=data.frame()
		outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
		outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	
		#输出train组风险文件
		riskScore=predict(multiCox,type="risk",newdata=train)           #利用train得到模型预测train样品风险
		coxGene=rownames(multiCoxSum$coefficients)
		coxGene=gsub("`","",coxGene)
		outCol=c("futime","fustat",coxGene)
		medianTrainRisk=median(riskScore)
		risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
		trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
		
		#输出test组风险文件
		riskScoreTest=predict(multiCox,type="risk",newdata=test)      #利用train得到模型预测test样品风险
		riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
		testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
		
		diff=survdiff(Surv(futime, fustat) ~risk,data = train)
		pValue=1-pchisq(diff$chisq,df=1)
		roc = survivalROC(Stime=train$futime, status=train$fustat, marker = riskScore, predict.time =1,  method="KM")
		
		diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
		pValueTest=1-pchisq(diffTest$chisq,df=1)
		rocTest = survivalROC(Stime=test$futime, status=test$fustat, marker = riskScoreTest, predict.time =1,  method="KM")
	
		if((pValue<0.05) & (roc$AUC>0.60) & (pValueTest<0.05) & (rocTest$AUC>0.60)){
		     #输出分组结果
			   write.table(trainOut,file="04.train.txt",sep="\t",quote=F,row.names=F)
			   write.table(testOut,file="04.test.txt",sep="\t",quote=F,row.names=F)
			   print('pvalue:')
			   print(pValue)
			   print('roc$AUC:')
			   print(roc$AUC)
			   print('pvalueTest:')
			   print(pValueTest)
			   print('rocTest$AUC:')
			   print(rocTest$AUC)
			   #输出单因素结果
			   write.table(outTab,file="05.uniCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(uniSigExp,file="05.uniSigExp.txt",sep="\t",row.names=F,quote=F)
			   #输出lasso回归结果
			   write.table(lassoSigExp,file="06.lassoSigExp.txt",sep="\t",row.names=F,quote=F)
			   pdf("06.lambda.pdf")
	       plot(fit, xvar = "lambda", label = TRUE)
	       dev.off()
	       pdf("06.cvfit.pdf")
	       plot(cvfit)
	       abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
	       dev.off()
	       #输出多因素结果
			   write.table(outMultiTab,file="07.multiCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(testRiskOut,file="riskTest.txt",sep="\t",quote=F,row.names=F)
			   write.table(trainRiskOut,file="riskTrain.txt",sep="\t",quote=F,row.names=F)
			   break
		}
}

#绘制森林图
options(forestplot_new_page = FALSE)
pdf(file="07.forest.pdf",width = 8,height = 5)
ggforest(multiCox,main = "Hazard ratio",cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

ibrary(survival)
library(survminer)
setwd("C:\\Users\\10982\\Desktop")       #设置工作目录

#定义生存曲线的函数
bioSurvival=function(inputFile=null,outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t")
	#比较高低风险组生存差异，得到显著性p值
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=TRUE,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
	print(surPlot)
	dev.off()
}
bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")

library(pheatmap)         #引用包
setwd("C:\\Users\\10982\\Desktop\\新建文件夹\\0.4\\model")      #设置工作目录
     #设置工作目录

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
	rt=rt[order(rt$riskScore),]      #按照riskScore对样品排序
		
	#绘制风险曲线
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=riskScoreFile, width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)", ylab="Risk score",
		 col=c(rep("green",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	dev.off()
		
	#绘制生存状态图
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="green"
	pdf(file=survStatFile, width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
		
	#绘制风险热图
	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(type=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=heatmapFile, width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("green",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}
#tarin组风险曲线
bioRiskPlot(inputFile="trainRisk.txt",riskScoreFile="train.riskScore.pdf",survStatFile="train.survStat.pdf",heatmapFile="train.heatmap.pdf")
#test组风险曲线
bioRiskPlot(inputFile="testRisk.txt",riskScoreFile="test.riskScore.pdf",survStatFile="test.survStat.pdf",heatmapFile="test.heatmap.pdf")

library(limma)
library(scatterplot3d)

setwd("C:\\Users\\10982\\Desktop")                        #???ù????

myPCA=function(input=null,output=null)
{
		#??????????????
		rt=read.table(input,sep="\t",header=T,check.names=F)
		rt=as.matrix(rt)
		rownames(rt)=rt[,1]
		exp=rt[,2:ncol(rt)]
		dimnames=list(rownames(exp),colnames(exp))
		data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		data=avereps(data)
		data=data[rowMeans(data)>0.5,]
		type=sapply(strsplit(colnames(data),"\\-"),"[",4)
		type=sapply(strsplit(type,""),"[",1)
		type=gsub("2","1",type)
		data=t(data[,type==0])
		rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
		#???risk???????
		risk=read.table("risk.txt",sep="\t",header=T,row.names=1)
		sameSample=intersect(rownames(data),rownames(risk))
		data=data[sameSample,]
		risk=risk[sameSample,]
		group=as.vector(risk[,"risk"])
		
		#PCA????
		data.class <- rownames(data)
		data.pca <- prcomp(data, scale. = TRUE)
		
		#?????
		color=ifelse(group=="low",3,2)
		pcaPredict=predict(data.pca)
		pdf(file=output,width=5.5,height=5)
		s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
		legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
		dev.off()
}

#??????л?????????????PCA???????04???symbol.txt???????
#myPCA(input="symbol.txt",output="allGene.PCA.pdf")
#????????????????????????PCA???????08???immuneGeneExp.txt???????
#myPCA(input="immuneGeneExp.txt"??output="immuneGene.PCA.pdf")
#???????????lncRNA??????????PCA???????09???immuneLncRNAexp.txt???????
#myPCA(input="immuneLncRNAexp.txt",output="immuneLncRNA.PCA.pdf")

#??????????lncRNA??????????PCA???????12???risk.txt???????
risk=read.table("testRisk.txt",sep="\t",header=T,row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])

data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
		
#?????
color=ifelse(group=="low",3,2)
pcaPredict=predict(data.pca)
pdf(file="riskTest.PCA.pdf",width=7,height=5)
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
dev.off()

###########################################################
risk=read.table("trainRisk.txt",sep="\t",header=T,row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])		
#PCA????
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
		
#?????
color=ifelse(group=="low",3,2)
pcaPredict=predict(data.pca)
pdf(file="RiskTrain.PCA.pdf",width=7,height=5)
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
dev.off()

#####################################
options(stringsAsFactors=F)
library(survival)          #引用包
setwd("D:\\biowolf\\GILnc")        #设置工作目录

#定义独立预后分析函数
indep=function(riskFile=null, cliFile=null, uniOutFile=null, multiOutFile=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
	
	#数据合并
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	
	#单因素独立预后分析
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
	
	
	#多因素独立预后分析
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
	rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
}

#独立预后分析
indep(riskFile="trainRisk.txt", cliFile="clinical.txt", uniOutFile="train.uniCox.txt", multiOutFile="train.multiCox.txt")
indep(riskFile="testRisk.txt", cliFile="clinical.txt", uniOutFile="test.uniCox.txt", multiOutFile="test.multiCox.txt")


library(survivalROC)
library(survival)
setwd("I:\\生信分析文章")      #设置工作目录
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取cox回归风险文件
rt$futime=rt$futime/365
rocCol=rainbow(3)
aucText=c()

model <- coxph(Surv(rt$futime,rt$fustat)~grade+stage+histological_type,data=rt)
beta <- summary(model)$coeff[,1]
rt$clinicalScore <- as.matrix(rt[,3:5])%*%beta

model <- coxph(Surv(rt$futime,rt$fustat)~grade+stage+histological_type+riskScore,data=rt)
beta <- summary(model)$coeff[,1]
rt$integration <- as.matrix(rt[,3:6])%*%beta

#绘制risk score的ROC曲线
pdf(file="multiROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")")
abline(0,1)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$clinicalScore, predict.time =5, method="KM")
aucText=c(aucText,paste0("clinical factor"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)

roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$integration, predict.time =5, method="KM")
aucText=c(aucText,paste0("clinical factor + risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[3],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

library(survivalROC)
setwd("I:\\生信分析文章")                     #设置工作目录
rt=read.table("indepInput.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取cox回归风险文件
rt$futime=rt$futime/365
rocCol=rainbow(ncol(rt)-2)
aucText=c()

#绘制risk score的ROC曲线
pdf(file="multiROC5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="False positive rate", ylab="True positive rate",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制其他临床性状的ROC曲线
j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
	roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =5, method="KM")
	j=j+1
	aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
	lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

library(rms)
setwd("C:\\Users\\10982\\Desktop")              #设置工作目录

#列线图绘制
riskFile="tcgaRisk.txt"
cliFile="tcgaClinical.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          #读取临床文件
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
#数据打包
dd <- datadist(rt)
options(datadist="dd")
#生成函数
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
#nomogram可视化
pdf(file=outFile,height=7.5,width=11)
plot(nom)
dev.off()

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("C:\\Users\\10982\\Desktop")           #设置工作目录
files=grep(".xls",dir(),value=T)                                         #获取目录下的所有xls文件
data = lapply(files,read.delim)                                          #读取每个文件
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                            #将文件后缀删掉

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
  theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
  theme(legend.background = element_blank()) + theme(legend.key = element_blank())
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "High risk<----------->Low risk", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

#将图形可视化，保存在"multipleGSEA.pdf"
pdf('multipleGSEA.pdf',      #输出图片的文件
     width=9,                #设置输出图片高度
     height=5)               #设置输出图片高度
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()

setwd("C:\\Users\\10982\\Desktop") 
library(pRRophetic)
library(ggplot2)
library(cowplot)
dat <- read.table("low-high.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
dat[1:3, 1:3]

ann <- read.table("type.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
head(ann)

table(ann$ImmClust)
GCP.drug <- read.table("drug.txt") 
GCP.drug <- GCP.drug$V1


jco <- c("#EABF00", "#2874C5", "red")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) 
  cat(drug," starts!\n") 
  

  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) 
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} 
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "ImmClust"=ifelse(ann$ImmClust == "low","Low-risk","High-risk"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("Low-risk","High-risk"),ordered = T) 
  
  # 绘图
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = ImmClust)) + 
    scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
    theme(legend.position="none") + # 倾斜字体
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p
  cat(drug," has been finished!\n") 
}


# 适合展示多种药物
p2 <- plot_grid(plotlist=plotp, ncol=6)
ggsave("boxplot of predicted IC50_multiple.pdf", width = 45, height =30)
###################计算p值
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "Low-risk"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "High-risk"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}
names(p) <- GCP.drug
print(p) 

write.table(p,"output_pvalue.txt", quote = F, sep = "\t")
