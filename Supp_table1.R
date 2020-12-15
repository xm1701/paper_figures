library(MatrixEQTL)
library(MASS)
################
n.sample<-matrix(c(500,500,10000,10000),2,2) 
n.sample.third<-c(2000)
n.all.snps<-10000
n.causal.snps<-n.all.snps*0.2
prop.m1.m2<-c(1)
###############
n.null.snps<-10000-n.causal.snps 
effect_mean<-c(0,0)
effect_var<-c(1)
effect_cor<-c(0,0.3,0.4)
effect_cor_part<-c(0.05,0.2,0.8)
###############
n.run<-500
###
CORR<-array(NA,c(n.run,2,3,3)) 
###
###############################################################
###############################################################
###############################################################
cat(" - alpha and beta, two set of summary statistics.\n")
set.seed(1234+2000+200*2) 
################################################################
############ 1.Generate Data1/2/3 SNPs##########################
################################################################
###########
(p.null<-n.null.snps)
(n.train<-n.sample[1,2])
x.train.null0<-matrix(NA,n.train,p.null)
subn.train<-n.train/5
############
(n.test<-n.sample[2,2])
x.test.null0<-matrix(NA,n.test,p.null)
###########NULL SNPs###########
for(i in 1:p.null){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###train data 
  x.train.null0[1:n.train,i]<-rbinom(n.train,2,mfa1)
  ###test data 
  x.test.null0[1:n.test,i]<-rbinom(n.test,2,mfa1)
}
###########CAUSAL SNPs###########
p.causal<-n.causal.snps
x.train.causal0<-matrix(NA,n.train,p.causal)
x.test.causal0<-matrix(NA,n.test,p.causal)
for(i in 1:p.causal){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###train data 
  x.train.causal0[1:n.train,i]<-rbinom(n.train,2,mfa1)
  ###testing data 
  x.test.causal0[1:n.test,i]<-rbinom(n.test,2,mfa1)
}
#######third SNP data###########
n.third<-n.sample.third
x.third0<-matrix(NA,n.third,n.all.snps)
for(i in 1:n.all.snps){
  mfa1<-runif(n=1, min = 0.05, max = 0.45)
  ###third data 
  x.third0[1:n.third,i]<-rbinom(n.third,2,mfa1)
}
##
x.train.null0[,which((apply(x.train.null0, 2, var))==0)]<-rbinom(n.train,2,0.5) 
x.train.causal0[,which((apply(x.train.causal0, 2, var))==0)]<-rbinom(n.train,2,0.5)
x.test.null0[,which((apply(x.test.null0, 2, var))==0)]<-rbinom(n.test,2,0.5)
x.test.causal0[,which((apply(x.test.causal0, 2, var))==0)]<-rbinom(n.test,2,0.5)
x.third0[,which((apply(x.third0, 2, var))==0)]<-rbinom(n.third,2,0.5)
##
x.train.null<-scale(x.train.null0)
x.test.null<-scale(x.test.null0)
x.train.causal<-scale(x.train.causal0)
x.test.causal<-scale(x.test.causal0)
x.third<-scale(x.third0)
##
x.train<-cbind(x.train.causal,x.train.null)
x.test<-cbind(x.test.causal,x.test.null)
###############################################################
############ 2.Generate Phenotypes ############################
###############################################################
for (h in 1:length(effect_cor)){
  for(j in 1: length(effect_cor_part)){
    for (iter in 1:n.run){if(iter %%10==0) {print(10000+iter)}
      set.seed(1234+iter+2000*h+200*l) 
      #####
      etrain.temp<-rnorm(n.train,mean=0,sd=sqrt(effect_var))
      etest.temp<-rnorm(n.test,mean=0,sd=sqrt(effect_var))
      #####
      Sigma.beta<- matrix(c(1,effect_cor[h],
                            effect_cor[h],1), ncol=2)   
      p.causal<-n.all.snps*effect_cor_part[j]
      x.train.causal<-x.train[,1:p.causal]
      x.test.causal<-x.test[,1:p.causal]
      coef<-mvrnorm(n=p.causal,mu=effect_mean,Sigma=effect_var*Sigma.beta)
      cor(coef)
      #####
      beta1<-as.matrix(coef[,1])
      beta2<-as.matrix(coef[,2])
      #####
      ###################################
      signal.train<-scale(x.train.causal%*%beta1)
      signal.test<-scale(x.test.causal%*%beta2)
      var(signal.train)
      ###################################
      y.train<-signal.train+etrain.temp
      #####
      y.test<-signal.test+etest.temp
      dim(y.train)
      dim(y.test)
      #####
      h2<-0.5
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
      ####Pull Matrix eQTL results
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
      ####Pull Matrix eQTL results
      pvalue.temp2<-m.temp2$all$eqtls$pvalue
      snps.temp02<-m.temp2$all$eqtls$snps
      snps.temp2<-as.numeric(substring(snps.temp02, 4))
      beta.temp2<-m.temp2$all$eqtls$beta
      eqtl.temp2<-data.frame(pvalue.temp2,snps.temp2,beta.temp2)
      ####################testing on the third############################
      ###all SNPs 
      test.pred1<-x.third[,eqtl.temp$snps.temp]%*%as.matrix(eqtl.temp$beta.temp)
      ###all SNPs 
      test.pred2<-x.third[,eqtl.temp2$snps.temp2]%*%as.matrix(eqtl.temp2$beta.temp2)
      ## 
      (CORR[iter,1,h,j]<-summary(lm(test.pred1~test.pred2))$coef[2,4])
      ####################testing on the third############################
      test.pred11<-x.test[,eqtl.temp$snps.temp]%*%as.matrix(eqtl.temp$beta.temp)
      ##
      (CORR[iter,2,h,j]<-summary(lm(y.test~test.pred11))$coef[2,4])
      ##
    }
  }
}
  
