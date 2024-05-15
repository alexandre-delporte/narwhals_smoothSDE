# HEADER --------------------------------------------
#
# Author:     Alexandre Delporte
# Copyright     Copyright 2024 - Alexandre Delporte
# Email:      alexandre.delporte@univ-grenoble-alpes.fr
#
# Date:     2024-05-05
#
# Script Name:    get_results.R
#
# Script Description: Compute biais,rmse and histograms from the results in the .csv files in the domains
# subdirectories
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(ggplot2)


domains=c("fjords")
times=c(12,24)
covs=c("DistanceShore")

for (domain in domains) {
  
  for (cov in covs) {
  
    for (time in times) {
  
      # Read csv files ----------
      pattern<- sprintf("result_%dh_%s_%s[1-9][0-9]*", time, domain,cov)
      files<- list.files(path=domain,pattern = pattern)
  
    
      #bind all results in one dataframe
      est_df=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
      
      for (file in files) {
        df=read.csv(file.path(domain,file))
        est_df=rbind(est_df,df)
      }
  
      # Bias and rmse for spline coefficients ------------
  
      est_df=est_df[-1,]
  
      N=length(files)
  
      biais=c()
      re_biais=c()
      rmse=c()
      true_values=c()
  
      #loop over the coefficient names
      for (name in unique(est_df$coeff_name)) {
    
        #all estimates for this coeff
        estimates=est_df[est_df$coeff_name==name,"estimate"]
        true_value=est_df[est_df$coeff_name==name,"true"][1]
    
        true_valuese=c(true_values,true_value)
    
        if (is.na(true_value)) {
          re_biais=c(re_biais,NA)
          biais=c(biais,NA)
          rmse=c(rmse,NA)
        }
        else {
      
          biais=c(biais,1/N*sum(estimates-true_value))
          rmse=c(rmse,sqrt(1/N*sum((estimates-true_value)^2)))
      
          if (true_value==0) {
            re_biais=c(re_biais,NA)
          }
          else {
            re_biais=c(re_biais,100*1/N*sum((estimates-true_value)/true_value))
          }
        }
      
        estimates=data.frame("estimate"=estimates)
        histo=ggplot(estimates,aes(x=estimate))+ geom_histogram(aes(y=..density..), colour="black", fill="white")+
          geom_density(alpha=.2, fill="#FF6666") +geom_vline(aes(xintercept=true_value),
                    color="blue", linetype="dashed", size=1)
      
        ggsave(filename=paste(paste("histo",name,domain,time,cov,sep="_"),".png",sep=""),plot=histo,width=10,height=5)
      }
  
      results=data.frame("coeff_name"=unique(est_df$coeff_name),"true"=true_values,
                                   "biais"=biais,"re_biais"=re_biais,"rmse"=rmse)
  
  
      p=ggplot()+geom_violin(data=est_df,aes(x=coeff_name,y=estimate,fill=coeff_name))+
        xlab(" ")+labs(fill = "Estimated coefficients")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        geom_point(data = df_DistanceShore, aes(x = coeff_name, y = true), color = "red", size = 1,shape=4)
  
  
      #round values to 2 digits
      round_numeric <- function(x,digits=3) {
        if (is.numeric(x)) {
          return(round(x, digits = digits))
        } else {
          return(x)
        }
      }
  
      results=as.data.frame(lapply(results, round_numeric))
  
      write.csv(results_DistanceShore, paste(paste("result_",domain,cov,time,sep="_"),".csv",sep=""), row.names=FALSE)
      }
  }
}