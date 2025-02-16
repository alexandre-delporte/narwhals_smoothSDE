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
library(dplyr)
library(plotly)
library(htmlwidgets)
library(here)


domains=c("fjords")

study_type="Ishore"

for (domain in domains) {
  
  hyperparams_files_list=list.files(path=domain,pattern = ".*.txt")
  
  for (hyperparams_file in hyperparams_files_list) {
    
    # Read the file and evaluate each line
    lines <- readLines(file.path(domain,hyperparams_file))
    for (line in lines) {
      eval(parse(text = line))
    }
      
    # get name of hyperparameters file without extension
    hyperparams_file_name=sapply(strsplit(hyperparams_file,"\\."),'[',1)
    
    
    # Read .out files and get percentage of trajectories that reached land
    out_files<- list.files(path=file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name)),
                           pattern ="\\.out$")
    
    # Initialize total trajectories count
    total_trajectories <- 0
    
    # Process each file
    for (file in out_files) {
      # Read the content of the file
      file_content <- readLines(file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name),file))
      
      # Find the line containing the percentage
      percentage_line <- grep("percent of the samples reached land", file_content, value = TRUE)
      
      if (length(percentage_line) > 0) {
        # Extract the percentage value using a regular expression
        percentage <- as.numeric(sub(".*([0-9.]+) percent of the samples reached land.*", "\\1", percentage_line))
        
        # Calculate the number of trajectories (12 samples per file)
        trajectories <- percentage / 100 * 12
        
        # Add to the total
        total_trajectories <- total_trajectories + trajectories
      }
    }
    
    # Print the total number of trajectories that reached land
    cat("Total percentage of trajectories that reached land:", total_trajectories/N_ID_HIGH/length(out_files)*100, "\n")
    
    # Read csv files ----------
    pattern<- sprintf("^estimates.*\\.csv$")
    result_files<- list.files(path=file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name)),pattern =pattern)
    

    #bind all results for each framework in a list of dataframes 
    all_results=list()
    
    #loop over output result files
    for (file in result_files) {
      
      # get framework type
      split<- strsplit(file, "_")[[1]]
      type=paste(split[2:5], collapse = "_")
      
      df=read.csv(file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name,sep=""),file))
      
      # remove rows for spline penalties and mu(fixed in the simulation study)
      
      df=df[!df$coeff_name %in% c("omega.te(theta,Ishore)","mu1.(Intercept)","mu2.(Intercept)"),]
      
      unstable=is.na(df[df$coeff_name=="tau.(Intercept)",3]) ||
        is.na(df[df$coeff_name=="nu.(Intercept)",3]) ||
        is.na(df[df$coeff_name=="tau.s(ID)",3]) ||
        is.na(df[df$coeff_name=="nu.s(ID)",3])
      
      if (unstable) {
        
        print(df[df$coeff_name=="tau.(Intercept)",2])
        next
      }
      
      if (type %in% names(all_results)) {
        all_results[[type]]<- rbind(all_results[[type]],df)
      }
      else {
        all_results[[type]]<- df
      }
    
    }
      
     
    #plot histograms for each parameter in each framework
    for (type in names(all_results)) {
      
      df=all_results[[type]]
      
      #get number of ID for this framework
      last_chain=sapply(strsplit(type,"_"),'[',4)
      N_ID=ifelse(last_chain=="lID",N_ID_LOW,N_ID_HIGH)
      
      # histogram for random effects in tau
      histo_tau_ID=ggplot(df[df$coeff_name %in% paste("tau.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        facet_wrap(~coeff_name)
      
      # histogram for random effects in nu
      histo_nu_ID=ggplot(df[df$coeff_name %in% paste("nu.s(ID).",1:N_ID,sep=""),],aes(x=estimate))+ 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        facet_wrap(~coeff_name)
      
      # histogram for standard deviation of random effects
      histo_sigma_re=ggplot(df[df$coeff_name %in% c("tau.s(ID)","nu.s(ID)"),],aes(x=estimate))+ 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        facet_wrap(~coeff_name)
      
      # histogram for intercepts tau and nu
      histo_intercepts=ggplot(df[df$coeff_name %in% c("tau.(Intercept)","nu.(Intercept)"),],aes(x=estimate))+ 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        facet_wrap(~coeff_name)
      
      # histogram for measurement error standard deviation
      histo_sigma_obs=ggplot(df[df$coeff_name == "log_sigma_obs",],aes(x=estimate))+ 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        facet_wrap(~coeff_name)
      
      

      
      ggsave(path=file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep="")),
                            filename=paste("histo_",type,"_tauID",".pdf",sep=""),
             plot=histo_tau_ID,width=10,height=5)
      ggsave(path=file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep="")),
                            filename=paste("histo_",type,"_nuID",".pdf",sep=""),
             plot=histo_nu_ID,width=10,height=5)
      ggsave(path=file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep="")),
             filename=paste("histo_",type,"_sigma_re",".pdf",sep=""),
             plot=histo_sigma_re,width=10,height=5)
      ggsave(path=file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep="")),
             filename=paste("histo_",type,"_intercepts",".pdf",sep=""),
             plot=histo_intercepts,width=10,height=5)
      ggsave(path=file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep="")),
      filename=paste("histo_",type,"_sigma_obs",".pdf",sep=""),
             plot=histo_sigma_obs,width=10,height=5)
      
    }
    
    
    # Mean, std error and CI -----------
    
    all_parameter_estimates=list()
    
    n_sim_list=list()
    
     #Loop over the frameworks 
    for (type in names(all_results)) {
      
      # results for single framework
      df=all_results[[type]]
      
      #keep only intercepts and standard deviations parameters
      df=df[df$coeff_name %in% c("tau.(Intercept)","nu.(Intercept)","tau.s(ID)","nu.s(ID)","log_sigma_obs"),]
      
      #name of coefficients
      coeff_names=unique(df$coeff_name)
      #link functions
      links=list(exp,exp,\(x) 1/sqrt(exp(x)),\(x) 1/sqrt(exp(x)),\(x) exp(x)*1000)
      names(links)=coeff_names
      #number of coeffs to estimate in single framework
      n_coeffs=length(coeff_names)
      # dataframe for estimates
      parameter_estimates_df=data.frame("mean"=rep(0,n_coeffs),
                                        "sd"=rep(0,n_coeffs),
                                        "low"=rep(0,n_coeffs),
                                        "up"=rep(0,n_coeffs))
      #loop over the coefficient names
      for (i in (1:n_coeffs)) {
      
        coeff_name=coeff_names[i]
        
        #all estimates for this coeff
        fn=links[[coeff_name]]
        estimates=fn(df[df$coeff_name==coeff_name,"estimate"])
        n_sim=length(estimates)
        
        if (!(type %in% names(n_sim_list))) {
          n_sim_list[[type]]=n_sim
        }
        
        mean_value=1/n_sim*sum(estimates)
        sd_value=sqrt(1/n_sim*sum((estimates-mean_value)^2))
        
        quantiles=quantile(estimates,c(0.025,0.975))
        
        parameter_estimates_df[i,]=c(mean_value,sd_value,quantiles)
    
      }
      
      coeff_names[5]="sigma_obs"
      parameter_estimates_df=cbind(coeff_names,parameter_estimates_df)
      
      
      all_parameter_estimates[[type]]<-parameter_estimates_df
    
    }
    
    print(n_sim_list)
     
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
      
      df=all_parameter_estimates[[type]]
      #save results in csv file
      output_file <- file.path(domain,paste("results_",study_type,"_",hyperparams_file_name,sep=""),
                               paste0("final_estimates_", type,".csv"))
      write.csv(df,output_file,
                row.names=FALSE)
    }
    
    
    # Plot the estimated smooth surfaces for omega
    
    # Read csv files ----------
    pattern<- sprintf("^surface.*\\.csv$")
    surface_files<- list.files(path=
                                 file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name)),
                               pattern =pattern)
    
    #a list of list to store all the estimated surfaces in the different scenarios
    all_surfaces=list()
    for (file in surface_files) {
      
      # get framework type
      split<- strsplit(file, "_")[[1]]
      type=paste(split[2:5], collapse = "_")
      
      df=read.csv(file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name,sep=""),file))

      
      
      if (type %in% names(all_surfaces)) {
        all_surfaces[[type]]<- c(all_surfaces[[type]],list(df))
      }
      else {
        all_surfaces[[type]]<- list(df)
      }
      
    }
    
    for (type in names(all_surfaces)) {
      
      
      estimated_surfaces=all_surfaces[[type]]
      
      # function to convert a dataframe to a matrix for surface plotting
      make_omega_matrix <- function(df) {
      unique_Dshore <- unique(df[df$Dshore>0.2 & df$Dshore<2,"Dshore"])
        unique_theta <- unique(df$theta)
        
        # Create a grid and fill matrix
        omega_matrix <- matrix(
          NA,
          nrow = length(unique_Dshore),
          ncol = length(unique_theta),
          dimnames = list(unique_Dshore, unique_theta)
        )
        
        for (i in seq_along(df$Dshore)) {
          dshore <- df$Dshore[i]
          theta <- df$theta[i]
          if (dshore>0.2 & dshore<2) {
            omega_matrix[as.character(dshore), as.character(theta)] <- df$omega[i]
          }
        }
        
        list(
          x = unique_Dshore,
          y = unique_theta,
          z = omega_matrix
        )
      }
      
      # Prepare data for each surface
      surfaces_matrices <- lapply(estimated_surfaces, make_omega_matrix)
      
      # Create a plotly object
      plot <- plot_ly()
      
      # Add each estimated surface to the plot
      for (i in seq_along(surfaces_matrices)) {
        surface <- surfaces_matrices[[i]]
        plot <- plot %>%
          add_surface(
            x = surface$x,
            y = surface$y,
            z = surface$z,
            opacity = 0.4,  
            showscale = FALSE,
            colorscale = list(c(0, "lightblue"), c(1, "lightblue"))
          )
      }
      
      fomega_cubic=function(cov_data,a,D0,D1,sigma_theta,sigma_D,b){
        Dshore=cov_data$DistanceShore
        theta=cov_data$theta
        
        a*theta*(theta-pi/2)*(theta+pi/2)*exp(-Dshore/D0)/Dshore+
          b*(exp(-1/2*(((theta+pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2))-
               exp(-1/2*(((theta-pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2)))
      }
      
      
      fomega=function(cov_data) {fomega_cubic(cov_data,A,D0,D1,SIGMA_THETA,SIGMA_D,B)} 
      
      # Smooth omega
      Dshore <- seq(from = 0.2, to = 2, length.out = 30)
      theta <- seq(from = -pi + 0.1, to = pi - 0.1, length.out = 30)
      grid <- as.data.frame(expand.grid(theta, Dshore))
      colnames(grid) <- c("theta", "DistanceShore")
      # Evaluate the true surface function for each grid point
      omega <-matrix(fomega(grid),30,30)
      
      # Add the true surface in red
      plot <- plot %>%
        add_surface(
          x = Dshore,
          y = theta,
          z = omega,
          opacity = 1,      
          showscale = FALSE,  
          colorscale = list(c(0, "red"), c(1, "red"))
        )
      
      # Customize layout
      plot <- plot %>%
        layout(
          title = NULL,  # Remove the title
          scene = list(
            xaxis = list(
              title = "Distance Shore",
              titlefont = list(size = 20),  # Increase axis label size
              tickfont = list(size = 16)    # Increase tick size
            ),
            yaxis = list(
              title = "Theta",
              titlefont = list(size = 20),
              tickfont = list(size = 16)
            ),
            zaxis = list(
              title = "Omega",
              titlefont = list(size = 20),
              tickfont = list(size = 16)
            ),
            aspectratio = list(x = 1, y = 1, z = 0.7)  # Adjust the aspect ratio
          )
        )
      

      # Display the plot
      setwd(file.path(getwd(),
                      file.path(domain,paste0("results_",study_type,"_",hyperparams_file_name,sep=""))))
      saveWidget(plot,paste0("final_surfaces_", type,".html"))
      setwd(file.path(here(),"R","simulation_study_CRCVM","estimation_cluster"))
      
      
    }
    
    
      

  
      
      
    }
  }

