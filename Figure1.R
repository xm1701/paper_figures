#***********************************************************************
#*********************Figure1*******************************************
#***********************************************************************
rm(list=ls())
################
ptm<- proc.time()
################
################
library(MatrixEQTL)
library(MASS)
################
n.sample<-matrix(c(10000,10000),1,2) #training data sample size 
n.sample.third<-c(2000)  #testing data sample size 
n.all.snps<-10000
n.causal.snps<-n.all.snps*0.2
prop.m1.m2<-c(1)
###############
n.null.snps<-10000-n.causal.snps #number of null SNPs
effect_mean<-c(0,0)
effect_var<-c(1)/n.all.snps
effect_cor<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)/n.all.snps
###############
n.run<-200 #number of replicates of the simulation
###
CORR<-array(NA,c(n.run,4,9)) 
###
cat(" - Start simulation.\n")
###############################################################
###############################################################
###############################################################
set.seed(1234+2000+200) 
################################################################
############ 1.Generate Data1/2/3 SNPs##########################
################################################################
###########
(p.null<-n.null.snps)
(n.train<-n.sample[1,1])
x.train.null0<-matrix(NA,n.train,p.null)
############
(n.test<-n.sample[1,2])
x.test.null0<-matrix(NA,n.test,p.null)
############
###########NULL SNPs###########
for(i in 1:p.null){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###train data 
  x.train.null0[1:(n.train),i]<-rbinom(n.train,2,mfa1)
  ###test data 
  x.test.null0[1:( n.test),i]<-rbinom(n.test,2,mfa1)
}
dim(x.train.null0)
###########CAUSAL SNPs###########
(p.causal<-n.causal.snps)
x.train.causal0<-matrix(NA,n.train,p.causal)
###########
x.test.causal0<-matrix(NA,n.test,p.causal)
#####
for(i in 1:p.causal){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###
  ###train data 
  x.train.causal0[1:( n.train),i]<-rbinom(n.train,2,mfa1)
  ###testing data 
  x.test.causal0[1:(n.test),i]<-rbinom(n.test,2,mfa1)
}
#######testing SNP data###########
(n.third<-n.sample.third)
x.third0<-matrix(NA,n.third,n.all.snps)
for(i in 1:n.all.snps){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###third data 
  x.third0[1:( n.third),i]<-rbinom(n.third,2,mfa1)
}
x.train.null0[,which((apply(x.train.null0, 2, var))==0)]<-rbinom(n.train,2,0.5) 
x.train.causal0[,which((apply(x.train.causal0, 2, var))==0)]<-rbinom(n.train,2,0.5)
x.test.null0[,which((apply(x.test.null0, 2, var))==0)]<-rbinom(n.test,2,0.5)
x.test.causal0[,which((apply(x.test.causal0, 2, var))==0)]<-rbinom(n.test,2,0.5)
x.third0[,which((apply(x.third0, 2, var))==0)]<-rbinom(n.third,2,0.5)
#######standardsize###########
##
x.train.null<-scale(x.train.null0)
x.test.null<-scale(x.test.null0)
##
x.train.causal<-scale(x.train.causal0)
x.test.causal<-scale(x.test.causal0)
##
x.third<-scale(x.third0)
##
dim(x.train.null);
dim(x.train.causal)
##
dim(x.test.null);
dim(x.test.causal)
##
dim(x.third)
x.train<-cbind(x.train.causal,x.train.null)
x.test<-cbind(x.test.causal,x.test.null)
###############################################################
############ 2.Generate Phenotypes ############################
###############################################################
###
iter<-1;h<-5; #an example for genetic correlation 0.5
###
for (h in 1:length(effect_cor)){
  for (iter in 1:n.run){if(iter %%10==0) {cat(paste(" - Simulation Iter: ",iter,"/",n.run,".",sep=""))}
    set.seed(1234+iter+2000*h+200) 
    #####
    Sigma.beta<- matrix(c(effect_var,effect_cor[h],
                          effect_cor[h],effect_var), ncol=2)   
    coef<-mvrnorm(n=p.causal,mu=effect_mean,Sigma=effect_var*Sigma.beta)
    cor(coef)
    #####
    beta1<-as.matrix(coef[,1])
    beta2<-as.matrix(coef[,2])
    #####heritability is set to one 
    y.train<-x.train.causal%*%beta1
    #####
    y.test<-x.test.causal%*%beta2
    ########################################################
    ############## 3. GWAS##################################
    ########################################################
    dim(x.train)
    ####################training data1############################
    ####Create three SlicedData objects for the analysis
    snps1 <- SlicedData$new(t(x.train));
    gene1 <- SlicedData$new(t(y.train));
    cvrt1 <- SlicedData$new();
    ####Produce no output files
    filename <- NULL; # tempfile()
    ####Call the main analysis function
    m.temp1 <- Matrix_eQTL_main(
      snps = snps1, 
      gene = gene1, 
      cvrt = cvrt1, 
      output_file_name = filename, 
      pvOutputThreshold = 1, 
      useModel = modelLINEAR, 
      errorCovariance = numeric(), 
      verbose = F,
      pvalue.hist = FALSE,
      min.pv.by.genesnp=F);
    ####Pull results
    pvalue.temp<-m.temp1$all$eqtls$pvalue
    snps.temp0<-m.temp1$all$eqtls$snps
    snps.temp<-as.numeric(substring(snps.temp0, 4))
    beta.temp<-m.temp1$all$eqtls$beta
    eqtl.temp<-data.frame(pvalue.temp,snps.temp,beta.temp)
    ######
    dim(x.test)
    ####################training data2############################
    ####Create three SlicedData objects for the analysis
    snps2 <- SlicedData$new(t(x.test));
    gene2 <- SlicedData$new(t(y.test));
    cvrt2 <- SlicedData$new();
    ####Produce no output files
    filename <- NULL; # tempfile()
    ####Call the main analysis function
    m.temp2 <- Matrix_eQTL_main(
      snps = snps2, 
      gene = gene2, 
      cvrt = cvrt2, 
      output_file_name = filename, 
      pvOutputThreshold = 1, 
      useModel = modelLINEAR, 
      errorCovariance = numeric(), 
      verbose = F,
      pvalue.hist = FALSE,
      min.pv.by.genesnp=F);
    ####Pull results
    pvalue.temp2<-m.temp2$all$eqtls$pvalue
    snps.temp02<-m.temp2$all$eqtls$snps
    snps.temp2<-as.numeric(substring(snps.temp02, 4))
    beta.temp2<-m.temp2$all$eqtls$beta
    eqtl.temp2<-data.frame(pvalue.temp2,snps.temp2,beta.temp2)
    ####################testing on the third dataset############################
    ###all SNPs 
    test.pred1<-x.third[,eqtl.temp$snps.temp]%*%as.matrix(eqtl.temp$beta.temp)
    ###all SNPs 
    test.pred2<-x.third[,eqtl.temp2$snps.temp2]%*%as.matrix(eqtl.temp2$beta.temp2)
    ##
    (CORR[iter,3,h]<-cor(test.pred1,test.pred2))
    (CORR[iter,4,h]<-cor(test.pred1,test.pred2)*sqrt((n.train+n.all.snps)*(n.test+n.all.snps)/(n.train*n.test)))
    ####################testing on the third dataset############################
    test.pred11<-x.test[,eqtl.temp$snps.temp]%*%as.matrix(eqtl.temp$beta.temp)
    #
    (CORR[iter,1,h]<-cor(test.pred11,y.test))
    (CORR[iter,2,h]<-cor(test.pred11,y.test)*sqrt((n.train+n.all.snps)/(n.train)))
  }
}


time.d1<-proc.time()-ptm
message(date(), ': Simulation finished. Time: about ',floor(time.d1[3]/60) , ' minute(s).')
################################################
cols<-c('cadetblue2',"deeppink1",'cadetblue2',"deeppink1")
ncase.pred<-c("A","B","C","D")
nylab<-c("Raw PRS Estimates","Corrected PRS Estimates","Raw PRS Estimates","Corrected PRS Estimates")
effect_cor2<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
##################
pdf("Figure-1.pdf",width=17,height=5)
########################
par(mar=c(5,5,5,5), xpd=F)
options(scipen=100000)
par(mfrow=c(1,4))
########################
for(ii in 1:4){
  #########################################
  #########################################
  boxplot((CORR[,ii,]), las = 1,
          col =cols[ii],
          cex.axis=1.5,cex.lab=2,cex.main=1.5,cex.sub=1.5,lwd=1.5,cex=2,
          xlab="True Genetic Correlation",ylab=paste(nylab[ii]),
          main=NULL, names =c(effect_cor2),
          ylim=c(0,1))
  grid(col = "lightgray", lty = "dotted",lwd=2,equilogs = TRUE)
  ##
  mtext(paste(ncase.pred[ii]), side = 3, line = 2,cex=1.5)
  ##
  abline(a=0,b=0.1,lty=3,lwd=4)
  #########################################
  #########################################
}


dev.off()

