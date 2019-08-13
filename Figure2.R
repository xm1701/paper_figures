#***********************************************************************
#*********************Figure2*******************************************
#***********************************************************************
rm(list=ls())
################
################
ptm<- proc.time()
################
################
library(MASS)
library(MatrixEQTL)
################
n.sample<-c(2000,10000) 
n.sample.test<-c(2000,2000)
n.causal.snps<-c(100,1000,5000,8000)
n.all.snps<-10000
n.null.snps<-n.all.snps-n.causal.snps #number of null SNPs
effect_mean<-c(0,0)
effect_var<-c(1)/n.all.snps
effect_cor<-c(0.8)/n.all.snps
##
n.run<-200 #number of replicates of the simulation
##
prs.cutoffs<-c(1,0.8,0.5,0.4,0.3,0.2,0.1,0.08,0.05,0.02,0.01,
               0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001)
n.prs<-length(prs.cutoffs)
###
CORR<-array(NA,c(n.run,2,4,n.prs)) 
###
###############################################################
###############################################################
###############################################################
iter<-1;l<-1;h<-4 # for example, 8000 causal SNPs, 2000 training samples  
###############################################################
###############################################################
###############################################################
for(l in 1:2){
  if(l==1){
    cat(" - Start simulation, small sampl size.\n")
  }
  ##
  if(l==2){
    cat(" - Start simulation, large sampl size.\n")
  }
  ##
  for(h in 1:length(n.null.snps)){
    set.seed(12345+2000*h+200*l) 
    ###############################################################
    ############ 1.Generate Data1/2 SNPs ##########################
    ###############################################################
    (p.null<-n.null.snps[h])
    (n.train<-n.sample[l])
    x.train.null0<-matrix(NA,n.train,p.null)
    (n.test<-n.sample.test[l])
    x.test.null0<-matrix(NA,n.test,p.null)
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
    (p.causal<-n.causal.snps[h])
    x.train.causal0<-matrix(NA,n.train,p.causal)
    ###########
    x.test.causal0<-matrix(NA,n.test,p.causal)
    #####
    for(i in 1:p.causal){
      mfa1<-runif(n=1, min = 0.05, max = 0.45)
      ###train data 
      x.train.causal0[1:( n.train),i]<-rbinom(n.train,2,mfa1)
      ###testing data 
      x.test.causal0[1:(n.test),i]<-rbinom(n.test,2,mfa1)
    }
    x.train.null0[,which((apply(x.train.null0, 2, var))==0)]<-rbinom(n.train,2,0.5) 
    x.train.causal0[,which((apply(x.train.causal0, 2, var))==0)]<-rbinom(n.train,2,0.5)
    x.test.null0[,which((apply(x.test.null0, 2, var))==0)]<-rbinom(n.test,2,0.5)
    x.test.causal0[,which((apply(x.test.causal0, 2, var))==0)]<-rbinom(n.test,2,0.5)
    #######standardsize###########
    ##
    x.train.null<-scale(x.train.null0)
    x.test.null<-scale(x.test.null0)
    ##
    x.train.causal<-scale(x.train.causal0)
    x.test.causal<-scale(x.test.causal0)
    ##
    dim(x.train.null);
    dim(x.train.causal)
    ##
    dim(x.test.null);
    dim(x.test.causal)
    x.train<-cbind(x.train.null,x.train.causal)
    x.test<-cbind(x.test.null,x.test.causal)
    ###############################################################
    ############ 2.Generate Phenotypes ############################
    ###############################################################
    for (iter in 1:n.run){if(iter %%10==0) {cat(paste(" - Simulation Iter: ",iter,"/",n.run,".",sep=""))}
      set.seed(12345+iter+2000*h+200*l) 
      #####
      Sigma.beta<- matrix(c(effect_var,effect_cor,
                            effect_cor,effect_var), ncol=2)   
      coef<-mvrnorm(n=p.causal,mu=effect_mean,Sigma=effect_var*Sigma.beta)
      cor(coef)
      #####
      beta1<-as.matrix(coef[,1])
      beta2<-as.matrix(coef[,2])
      #####
      y.train<-x.train.causal%*%beta1
      #####
      y.test<-x.test.causal%*%beta2
      ########################################################
      ############## 3. GWAS##################################
      ########################################################
      dim(x.train)
      ####################training############################
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
      ####################testing############################
      ###with thresholding 
      corr0.temp<-rep(NA,n.prs) #with y
      for(jj in 1:n.prs){
        ####
        select.temp<-eqtl.temp[eqtl.temp[,1]<prs.cutoffs[jj],2]
        select.beta.temp<-eqtl.temp[eqtl.temp[,1]<prs.cutoffs[jj],3]
        ####all selected 
        if(length(select.temp)<=1){
          corr0.temp[jj]<-NA
        }
        if(length(select.temp)>1){
          test.pred1<-x.test[,select.temp]%*%as.matrix(select.beta.temp)
          corr0.temp[jj]<-cor(test.pred1,y.test)
        }
      }
      ###
      (CORR[iter,l,h,]<-corr0.temp)
    }
    ##
  }
  ##########################
  ##########################
}

###############################################################
time.d1<-proc.time()-ptm
message(date(), ': Simulation finished. Time: about ',floor(time.d1[3]/60) , ' minute(s).')
###############################################################
###############################################################
################################################
ncase.pred<-c("A","B","C","D")
###########################boxplots####################################
cols<-c('orange','yellow','green','red','sienna','palevioletred1','royalblue2','darkturquoise',
        'darkviolet','brown','midnightblue','dimgrey','darkorange1',"forestgreen","blueviolet")
xtitles<-c("Cutoffs")
prs.cutoffs<-c(1,0.8,0.5,0.4,0.3,0.2,0.1,0.08,0.05,0.02,0.01,
               0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001)
name.pred<-list()
name.pred[[1]]<-(prs.cutoffs)
name.pred[[2]]<-c("A","B","C","D")
name.pred[[3]]<-c(0.01,0.1,0.5,0.8)
##################
pdf("Figure-2.pdf",width=15,height=8)
########################
par(mar=c(5.5,4.5,5,4), xpd=F)
par(mfrow=c(2,4))
#options(scipen=1000)
########################
for(i in c(2,1)){
  for(k in 1:4){
    #########################################
    #########################################
    boxplot((CORR[,i,k,]), las = 2,
            col =cols,lwd=1.5,cex=1.5,cex.lab=1.5,cex.main=1.5,
            names=prs.cutoffs,
            main=NULL,
            ylim=c(0,1),
            xlab=xtitles,
            ylab="Estimated Genetic Correlation")
    ##
    mtext(paste("m/p= ",name.pred[[3]][k],sep = ""), side = 3, line = 0)
    ##
    if(i==1){
      abline(h=sqrt(2000/(2000+10000))*0.8,col='forestgreen',lty=3,lwd=3)
      mtext(paste("n=",2000,", p=",10000,sep=""), side = 3, line = 1)
    }
    if(i==2){
      abline(h=sqrt(10000/(10000+10000))*0.8,col='forestgreen',lty=3,lwd=3)
      mtext(paste("n=",10000,", p=",10000,sep=""), side = 3, line = 1)
    }
    ##
    grid(col = "lightgray", lty = "dotted",lwd=3,equilogs = TRUE)
    abline(h=0.8,col="deeppink1",lty=3,lwd=3)
    #########################################
    #########################################
  }
}



dev.off()


