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
# Script Description: 
# Compute mean, biais,rmse and histograms from the results in the .csv files in the domains
# subdirectories
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(ggplot2)
library(tidyverse)


domains=c("fjords")
list_TMAX=c(12)
covs=c("DistanceShore")
list_N_ID=c(12)
list_DMIN=c(4)

for (domain in domains) {
  
  for (cov in covs) {
  
    for (TMAX in list_TMAX) {
      
      for (N_ID in list_N_ID) {
        
        for (DMIN in list_DMIN) {
  
          # Read csv files ----------
          pattern<- sprintf("result_%s_%dh_%dID_%dkm_%s[1-9][0-9]*",domain,TMAX,N_ID,DMIN,cov)
          files<- list.files(path=domain,pattern = pattern)
      

          #bind all results in one dataframe
          all_df=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
      
          for (file in files) {
            df=read.csv(file.path(domain,file))
            all_df=rbind(all_df,df)
          }
          
          # remove irrelevant rows
          all_df=all_df[-1,]
          all_df=subset(all_df,!(coeff_name %in% c("mu1.(Intercept)","mu2.(Intercept)","omega.te(theta,DistanceShore)")))
          
          #number of simulations and coeffs
          N=length(files)
          n_coeff=length(unique(all_df$coeff_name))
          
          # true values
          true_values=all_df[1:n_coeff,"true"]
          names(true_values)=unique(all_df$coeff_name)
  
          
          #plot histograms
          histo_tau_ID=ggplot(all_df[all_df$coeff_name %in% paste("tau.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
            geom_histogram(aes(y=..density..), colour="black", fill="white")+
            geom_density(alpha=.2, fill="#FF6666") +
            geom_vline(aes(xintercept=true),color="blue", linetype="dashed", size=0.5)+
            facet_wrap(~coeff_name)
          
          histo_nu_ID=ggplot(all_df[all_df$coeff_name %in% paste("nu.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
            geom_histogram(aes(y=..density..), colour="black", fill="white")+
            geom_density(alpha=.2, fill="#FF6666") +
            geom_vline(aes(xintercept=true),color="blue", linetype="dashed", size=0.5)+
            facet_wrap(~coeff_name)
          
          histo_omega_te=ggplot(all_df[all_df$coeff_name %in% paste(paste("omega.te(theta,",cov,").",sep=""),1:8,sep=""),],aes(x=estimate))+ 
            geom_histogram(aes(y=..density..), colour="black", fill="white")+
            geom_density(alpha=.2, fill="#FF6666") +
            geom_vline(aes(xintercept=true),color="blue", linetype="dashed", size=0.5)+
            facet_wrap(~coeff_name)
          
          histo_sigma=ggplot(all_df[all_df$coeff_name %in% c("tau.s(ID)","nu.s(ID)"),],aes(x=estimate))+ 
            geom_histogram(aes(y=..density..), colour="black", fill="white")+
            geom_density(alpha=.2, fill="#FF6666") +
            geom_vline(aes(xintercept=true),color="blue", linetype="dashed", size=0.5)+
            facet_wrap(~coeff_name)
          
          histo_intercepts=ggplot(all_df[all_df$coeff_name %in% c("tau.(Intercept)","nu.(Intercept)","omega.(Intercept)"),],aes(x=estimate))+ 
            geom_histogram(aes(y=..density..), colour="black", fill="white")+
            geom_density(alpha=.2, fill="#FF6666") +
            geom_vline(aes(xintercept=true),color="blue", linetype="dashed", size=0.5)+
            facet_wrap(~coeff_name)
        
          ggsave(filename=paste(paste("histo","tauID",domain,TMAX,cov,sep="_"),".png",sep=""),
                 plot=histo_tau_ID,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo","nuID",domain,TMAX,cov,sep="_"),".png",sep=""),
                 plot=histo_nu_ID,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo","omega.te",domain,TMAX,cov,sep="_"),".png",sep=""),
                 plot=histo_omega_te,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo","sigma",domain,TMAX,cov,sep="_"),".png",sep=""),
                 plot=histo_sigma,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo","intercepts",domain,TMAX,cov,sep="_"),".png",sep=""),
                 plot=histo_intercepts,path=file.path(domain),width=10,height=5)
  
          # Bias and rmse for spline coefficients ------------
  
          biais=c()
          re_biais=c()
          rmse=c()
          mean=c()
          sd=c()
      
          #loop over the coefficient names
          for (name in unique(all_df$coeff_name)) {
        
            #all estimates for this coeff
            estimates=all_df[all_df$coeff_name==name,"estimate"]
           
            #true value
            true_value=true_values[name]
         
            biais=c(biais,1/N*sum(estimates-true_value))
            rmse=c(rmse,sqrt(1/N*sum((estimates-true_value)^2)))
            mean_value=1/N*sum(estimates)
            mean=c(mean,mean_value)
            sd=c(sd,sqrt(1/N*sum((estimates-mean_value)^2)))
          
            if (true_value==0) {
              re_biais=c(re_biais,NA)
            }
            else {
                re_biais=c(re_biais,100*1/N*sum((estimates-true_value)/true_value))
            }
            }

          
          results=data.frame("coeff_name"=unique(all_df$coeff_name),"true"=true_values,"mean"=mean,
                             "sd"=sd,"biais"=biais,"re_biais"=re_biais,"rmse"=rmse)
      
      
          p=ggplot()+geom_violin(data=all_df,aes(x=coeff_name,y=estimate,fill=coeff_name))+
            xlab(" ")+labs(fill = "Estimated coefficients")+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            geom_point(data = all_df, aes(x = coeff_name, y = true), color = "red", size = 1,shape=4)
      
      
          #round values to 2 digits
          round_numeric <- function(x,digits=3) {
            if (is.numeric(x)) {
              return(round(x, digits = digits))
            } else {
              return(x)
            }
          }
      
          results=as.data.frame(lapply(results, round_numeric))
      
          write.csv(results, paste(paste("result_",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".csv",sep=""), row.names=FALSE)
        }
      }
    }
  }
}
