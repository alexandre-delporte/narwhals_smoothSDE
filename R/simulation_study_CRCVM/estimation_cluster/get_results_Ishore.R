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
    
     #Loop over the frameworks 
    for (type in names(all_results)) {
      
      # results for single framework
      df=all_results[[type]]
      
      #keep only intercepts and standard deviations parameters
      df=df[df$coeff_name %in% c("tau.(Intercept)","nu.(Intercept)","tau.s(ID)","nu.s(ID)","log_sigma_obs"),]
      
      #name of coefficients
      coeff_names=unique(df$coeff_name)
      #link functions
      links=list(exp,exp,identity,identity,\(x) exp(x)*1000)
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
        
        mean_value=1/n_sim*sum(estimates)
        sd_value=sqrt(1/n_sim*sum((estimates-mean_value)^2))
        
        quantiles=quantile(estimates,c(0.025,0.975))
        
        parameter_estimates_df[i,]=c(mean_value,sd_value,quantiles)
    
      }
      
      coeff_names[5]="sigma_obs"
      parameter_estimates_df=cbind(coeff_names,parameter_estimates_df)
      
      
      all_parameter_estimates[[type]]<-parameter_estimates_df
    
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
      
      print(head(df))
      
      
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
        unique_Dshore <- unique(df$Dshore)
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
          omega_matrix[as.character(dshore), as.character(theta)] <- df$omega[i]
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
            opacity = 0.1,       # Make surfaces semi-transparent
            showscale = FALSE    # Hide color scale for individual surfaces
          )
      }
      
      # Customize layout
      plot <- plot %>%
        layout(
          title = "3D Plot of Estimated Surfaces and Mean Surface",
          scene = list(
            xaxis = list(title = "Distance Shore"),
            yaxis = list(title = "Theta"),
            zaxis = list(title = "Omega")
          )
        )
      
      # Display the plot
      saveWidget(plot,paste0("final_surfaces_", type,".html"))
      
      
    }
    
    
      

  
      
      
    }
  }

