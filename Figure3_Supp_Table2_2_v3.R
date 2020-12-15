
####################################################################
##############################PING##################################
####################################################################
####read in the data of seven ROI volumes:thalamus proper, caudate, putamen, pallidum, hippocampus, accumbens area, and the total brain volume (TBV).
pheno7<-read.table("ping_roi_volume_seven.phen")
####read in the coovariates data (age and gender)
covs<-read.csv("PING_Covariates.csv",header =T)
##########################
##########################
summary.reaction02<-array(NA,c(7,2,11))
for(j in 1:7){
  f1_score<-read.table(paste0("Reaction_PING_c1.profile"),header =T)
  dim(f1_score)
  pheno.true<-rep(NA,length(f1_score$FID))
  cov.true<-matrix(NA,length(f1_score$FID),2)
  ###############
  for (i in 1:length(f1_score$FID)){
    temp.id<-which(as.character(pheno$ID)==as.character(f1_score$IID[i]))
    if(length(temp.id)>0){
      pheno.true[i]<-pheno7[temp.id[1],j]
    }
    temp.id2<-which(as.character(covs$ID)==as.character(f1_score$IID[i]))
    if(length(temp.id2)>0){
      cov.true[i,1]<-covs$age[temp.id2]
      cov.true[i,2]<-covs$gender[temp.id2]
    }
  }
  ###############
  for(i in 1:11){
    ######reaction#######
    f1_score<-read.table(paste0("Reaction_PING_c",i,".profile"),header =T)
    if (class(f1_score) == "try-error"){
      next;
    }
    dim(f1_score)
    l1<-lm(pheno.true~f1_score$SCORE+cov.true)
    l0<-lm(pheno.true~cov.true)
    (summary.reaction02[j,1,i]<-summary(l1)$coef[2,4])# the p-value 
    (summary.reaction02[j,2,i]<-summary(l1)$r.squared-summary(l0)$r.squared)# the partial R^2
    ###############
  }
}

################################################
##############Generate Figure 3#################
################################################
cols<-c('hotpink','royalblue2','darkturquoise',
        'brown',"red",'dimgrey','darkorange1',"forestgreen","blueviolet",'orange')
pdf("Figure3.pdf",width=10,height=8)
########################
par(mar=c(4,4,4,4), xpd=F)
par(mfrow=c(1,1))
for(i in 1:7){
  if(i==1){
    plot(y=summary.reaction02[1,2,]*100,x=c(1:11),ylim = c(0,4.5),type="b",col=cols[i],lwd=4,lty=1,
         ylab = "",xlab = "",cex=1, xaxt = "n")
    axis(1, at=1:11, labels=c(1,0.8,0.5,0.4,0.3,0.2,0.1,0.08,0.05,0.02,0.01))
    grid(col = "lightgray", lty = "dotted",lwd=3,equilogs = TRUE)
    mtext(paste("Cutoffs"), side = 1, line = 2.5,cex=1.5)
    mtext(paste("Raw Partial R^2 (X100%)"), side = 2, line = 2.5,cex=1.5)
    #########
    legend("topleft", legend=c("thalamus proper","caudate","putamen","pallidum", "hippocampus","accumbens area", "total brain volume"),
           lty=c(1),
           lwd=4,
           col=c(cols[1:7]), 
           inset=-0,cex = 1.4,bty ="n")
    #########
  }
  if(i>1){
    lines(summary.reaction.ping2[i,2,]*100,col=cols[i],lwd=4,type="b")
  }
}
dev.off()



