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
hyperparams_files_list=list.files(pattern = ".*.txt")

for (domain in domains) {
  
  for (hyperparams_file in hyperparams_files_list) {
          
          # get name of hyperparameters file without extension
          hyperparams_file_name=sapply(strsplit(hyperparams_file,"\\."),'[',1)
  
          # Read csv files ----------
          pattern<- sprintf("^result_%s.*%s\\.csv$",domain,hyperparams_file_name)
          result_files<- list.files(path=domain,pattern =pattern)
      

          #bind all results for each framework in a list of dataframes 
          all_results=list()
          
          #loop over output result files
          for (file in result_files) {
            
            # get framework type
            split<- strsplit(file, "_")[[1]]
            type=paste(split[3:6], collapse = "_")
            
            df=read.csv(file.path(domain,file))
            
            if (type %in% names(all_results)) {
              all_results[[type]]<- rbind(all_results[[type]],df)
            }
            else {
              all_results[[type]]<- df
            }
          
          }
          
         
          
          #plot histograms for each parameter in each framework
          for (type in names(all_results)) {
            
            
            # histogram for random effects in tau
            histo_tau_ID=ggplot(all_results[[type]][all_results$coeff_name %in% paste("tau.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.2, fill="#FF6666") +
              facet_wrap(~coeff_name)
            
            # histogram for random effects in nu
            histo_nu_ID=ggplot(all_results[all_results$coeff_name %in% paste("nu.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.2, fill="#FF6666") +
              facet_wrap(~coeff_name)
            
            # histogram for standard deviation of random effects
            histo_sigma_re=ggplot(all_results[all_results$coeff_name %in% c("tau.s(ID)","nu.s(ID)"),],aes(x=estimate))+ 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.2, fill="#FF6666") +
              facet_wrap(~coeff_name)
            
            # histogram for intercepts tau and nu
            histo_intercepts=ggplot(all_results[all_results$coeff_name %in% c("tau.(Intercept)","nu.(Intercept)"),],aes(x=estimate))+ 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.2, fill="#FF6666") +
              facet_wrap(~coeff_name)
            
            # histogram for measurement error standard deviation
            histo_sigma_obs=ggplot(all_results[all_results$coeff_name == "log_sigma_obs",],aes(x=estimate))+ 
              geom_histogram(aes(y=..density..), colour="black", fill="white")+
              geom_density(alpha=.2, fill="#FF6666") +
              facet_wrap(~coeff_name)
            
            

            
            ggsave(filename=paste("histo_",domain,"_",type,"_tauID",".png",sep=""),
                   plot=histo_tau_ID,path=file.path(domain),width=10,height=5)
            ggsave(filename=paste("histo_",domain,"_",type,"_nuID",".png",sep=""),
                   plot=histo_nu_ID,path=file.path(domain),width=10,height=5)
            ggsave(filename=paste("histo_",domain,"_",type,"_sigma_re",".png",sep=""),
                   plot=histo_sigma_re,path=file.path(domain),width=10,height=5)
            ggsave(filename=paste("histo_",domain,"_",type,"_sigma_obs",".png",sep=""),
                   plot=histo_intercepts,path=file.path(domain),width=10,height=5)
            
          }
          
  
          # Mean, std error and CI -----------
          
          all_parameter_estimates=list()
          
           #Loop over the frameworks 
          for (type in names(all_results)) {
            
            # results for single framework
            df=all_results[[type]]
            #number of coeffs to estimate in single framework
            n_coeffs=length(unique(df$coeff_name))
            #number of simulations
            n_sim=length(df$coeff_name)/n_coeffs
            # dataframe for estimates
            parameter_estimates_df=data.frame("coeff_name"=unique(all_results[[type]])$coeff_name,
                                              "mean"=rep(0,n_coeffs),
                                              "sd"=rep(0,n_coeffs),
                                              "low"=rep(0,n_coeffs),
                                              "up"=rep(0,n_coeffs))
            
            #loop over the coefficient names
            for (coeff_name in unique(all_results[[type]]$coeff_name)) {
            
            #all estimates for this coeff
            estimates=df[df$coeff_name==coeff_name,"estimate"]
           
            mean_value=1/n_sim*sum(estimates)
            sd_value=sqrt(1/n_sim*sum((estimates-mean_value)^2))
            
            quantile=quantile(estimates,c(0.025,0.975))
            
            parameter_estimates_df[coeff_name,]=c(mean_value,sd_value,quantiles)
          
            }
            
          all_parameter_estimate[[type]]<-parameter_estimates_df
          
          }
         
          #round values to 2 digits
          round_numeric <- function(x,digits=3) {
            if (is.numeric(x)) {
              return(round(x, digits = digits))
            } else {
              return(x)
            }
          }
      
          all_parameter_estimates=lapply(all_parameter_estimates, function(df) {as.data.frame(lapply(df, round_numeric))})
          
          for (type in names(all_parameter_estimates)) {
            
            write.csv(path=file.path(domain),all_parameter_estimates[[type]], 
                      paste("result_",domain,"_",type,"_",hyperparams_file_name,".csv",sep=""),
                      row.names=FALSE)
          
          
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
