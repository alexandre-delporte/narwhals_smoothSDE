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
# Compute mean, std error, biais,rmse and histograms from the results in the .csv files in the domains
# subdirectories
#
#
# SETUP ------------------------------------
cat("\014")                 # Clears the console
rm(list = ls())             # Remove all variables of the work space

library(ggplot2)
library(tidyverse)


domains=c("circ","rect","fjords")
list_TMAX=c(12,60)
covs=c("DistanceShore")
list_N_ID=c(6,12)
list_DMIN=c(1)
SP_DF=c(3,3)
D_low=1
D_up=5

for (domain in domains) {
  
  for (cov in covs) {
  
    for (TMAX in list_TMAX) {
      
      print(TMAX)
      
      for (N_ID in list_N_ID) {
        
        for (DMIN in list_DMIN) {
  
          # Read csv files ----------
          pattern<- sprintf("result_%s_%dh_%dID_%dkm_%s[1-9][0-9]*",domain,TMAX,N_ID,DMIN,cov)
          files<- list.files(path=domain,pattern = pattern)
      

          #bind all results in one dataframe 
          all_df=data.frame("coeff_name"=NA,"estimate"=NA,"true"=NA)
          
          #loop over output result files
          for (file in files) {
            
            df=read.csv(file.path(domain,file))
            n_row=length(df$true)
            value=df[n_row-2,2]
            if (value<0.75) {
              all_df=rbind(all_df,df)
            }
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
    
          ggsave(filename=paste(paste("histo_tauID",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=histo_tau_ID,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo_nuID",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=histo_nu_ID,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo_omega.te",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=histo_omega_te,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo_sigma",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=histo_sigma,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("histo_intercepts",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=histo_intercepts,path=file.path(domain),width=10,height=5)
  
          # Mean, std error, bias and rmse ------------
  
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

          #create the dataframe with results
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
          
          
          #get the set up
          
          # Defintion of smooth parameter omega ----------
          
          fomega=function(cov_data,D0=0.5,omega0=60*pi,lambda=1.5,kappa=0.2) {
            Dshore=cov_data$DistanceShore
            theta=cov_data$theta
            if (is.null(Dshore)){
              Dshore=1/cov_data$ExpShore
            }
            coeff=exp(-kappa*(Dshore/D0)^2)
            omega=omega0/2*(tanh(lambda*(theta+pi/2-atanh(-0.9)/lambda))+tanh(lambda*(theta-(pi/2-atanh(-0.9)/lambda))))*coeff
            return(omega)
          }
          
          
          
          
          
          # Approximation of smooth omega with tensor splines -----------------------
          
          n <- 100000                           #number of observations
          D_low=1
          D_up=5
          
          theta <- runif(n,-pi,pi)            #sample theta
          DistanceShore <- runif(n,D_low,D_up)     #sample DistanceShore
          
          samples=data.frame(theta=theta,DistanceShore=DistanceShore)  
          
          #define grid of values
          theta_v <- seq(-pi,pi,length=30)
          Dshore_v<- seq(0.1,7,length=30)
          pr <- data.frame(theta=rep(theta_v,30),DistanceShore=rep(Dshore_v,rep(30,30)))
          
          # true values of the function over this grid
          truth <- matrix(fomega(pr),30,30)
          
          #points on the surface perturbed by gaussian noise
          f <- fomega(samples)
          y <- f+0.1*rnorm(n)
          
          #plot true function
          persp(theta_v,Dshore_v,truth);
          title("truth")
          
          
          #fit with bivariate splines te fo DistanceShore
          knots_DistanceShore=list(theta=seq(-pi,pi,len=SP_DF[1]),DistanceShore=seq(D_low,D_up,len=SP_DF[2]))
          m1 <- gam(y~te(theta,DistanceShore,k=SP_DF,bs="cs"),knots=knots_DistanceShore)
          
          #visualize
          vis.gam(m1);title("tensor product")
          
          df=data.frame("x"=theta,"z"=DistanceShore,"y"=y)
          fig <- plot_ly(df,type="scatter3d", x = ~x, y = ~z, z = ~y,mode="markers",marker=list(size=4))
          fig
          
          
          
          ExpShore=1/DistanceShore
          
          #fit with bivariate splines te of ExpShore
          knots_ExpShore=list(theta=seq(-pi,pi,len=SP_DF[1]),ExpShore=seq(1/D_up,1/D_low,len=SP_DF[2]))
          m2 <- gam(y~te(theta,ExpShore,k=SP_DF,bs="cs"),knots=knots_ExpShore)
          
          
          #visualize
          vis.gam(m2);title("tensor product")
          
          
          #get splines coefficients
          
          sp_coeff_Dshore=m1$coefficients
          sp_coeff_ExpShore=m2$coefficients
          
          
          #define spline smooth parameter omega
          fomega_splines=function(cov_data) {
            
            if (is.null(cov_data$ExpShore)){
              omega=predict(m1,newdata=cov_data)
            }
            else if (is.null(cov_data$DistanceShore)){
              omega=predict(m2,newdata=cov_data)
            }
            
            return(as.numeric(omega))
          }
          
          #k-th basis spline function for omega
          fomega_spline_basis=function(cov_data,k=5) {
            
            if (is.null(cov_data$ExpShore)){
              Xp=predict(m1, newdata = cov_data, type = "lpmatrix")
              
            }
            else if (is.null(cov_data$DistanceShore)){
              Xp=predict(m2, newdata = cov_data, type = "lpmatrix")
            }
            coeffs=rep(0,ncol(Xp))
            coeffs[k]=1
            return(as.numeric(Xp%*%coeffs))
          }
          
        
          angle1=pi/2
          angle2=pi
          # Generate the x-values
          x_vals <- seq(from=1.2, to=3.8, length.out=100)
          
          # True functions
          omega_away <- fomega_splines(cov_data=data.frame("theta"=rep(angle1, 100), "DistanceShore"=x_vals))
          omega_toward <- fomega_splines(cov_data=data.frame("theta"=rep(angle2, 100), "DistanceShore"=x_vals))
          
          # Create initial plots with true functions
          plot_omega_away <- ggplot() +
            geom_line(aes(x=x_vals, y=omega_away), color="red") +
            ylab(expression(omega)) +
            xlab("Distance to shore")
          
          plot_omega_toward <- ggplot() +
            geom_line(aes(x=x_vals, y=omega_toward), color="red") +
            ylab(expression(omega)) +
            xlab("Distance to shore")
          
          # Loop over the output files
          for (file in files) {
            df <- read.csv(file.path(domain, file))
            omega_coeffs <- as.numeric(df[df$coeff_name %in% paste("omega.te(theta,DistanceShore).", 1:(SP_DF[1]*SP_DF[2]-1), sep=""), "estimate"])
            
            fomega_estimate <- function(cov_data) {
              sum <- 0
              for (k in 1:(SP_DF[1]*SP_DF[2]-1)) {
                sum <- sum + omega_coeffs[k] * fomega_spline_basis(cov_data, k)
              }
              return(sum)
            }
            
            # Compute estimates
            est_omega_away <- fomega_estimate(cov_data=data.frame("theta"=rep(angle1, 100), "DistanceShore"=x_vals))
            est_omega_toward <- fomega_estimate(cov_data=data.frame("theta"=rep(angle2, 100), "DistanceShore"=x_vals))
            
            # Convert to data frames
            df_away <- data.frame(x=x_vals, y=est_omega_away)
            df_toward <- data.frame(x=x_vals, y=est_omega_toward)
            
            # Add estimates to the plots
            plot_omega_away <- plot_omega_away +
              geom_line(data=df_away, aes(x=x, y=y), color="lightblue")
            
            plot_omega_toward <- plot_omega_toward +
              geom_line(data=df_toward, aes(x=x, y=y), color="lightblue")
          }
          
          
          # Print plots
          print(TMAX)
          print(N_ID)
          print(paste(paste("omega_est_away_",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""))
          ggsave(filename=paste(paste("omega_est_away_",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=plot_omega_away,path=file.path(domain),width=10,height=5)
          ggsave(filename=paste(paste("omega_est_toward_",domain,"_",TMAX,"h_",N_ID,"ID_",DMIN,"km_",cov,sep=""),".png",sep=""),
                 plot=plot_omega_toward,path=file.path(domain),width=10,height=5)
        }
      }
    }
  }
}
