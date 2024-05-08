

#SEED FOR REPRODUCTIBILITY
set.seed(42)

library(smoothSDE)
library(ggplot2)
library(dplyr)
library(plotly)
library(htmlwidgets)




################################### OUTPUT FILE TO KEEP RECORD OF THE SCRIPT EXECUTION ####################
filename=paste("response_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

########################################    GET NARWHAL DATA   #############################################


# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")  


# DATA BEFORE EXPOSURE

dataBE1=read.csv(file.path(narwhal_data_path,"DataBE1.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)

dataBE2=read.csv(file.path(narwhal_data_path,"DataBE2.csv"), header = TRUE,dec = ".")



cat("Extracting trajectories before exposure and 12h after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


# DATA AFTER EXPOSURE

dataAE=read.csv(file.path(narwhal_data_path,"DataAE.csv"), header = TRUE,dec = ".")


cat("Extracting trajectories after exposure.. ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataAE[,1]),"positions measured after exposure \n"),file=filename,sep="\n",append=TRUE)


#######################  CIR  MODEL #######################

#define model
formulas <- list(mu=~1,
                 beta =~1,
                 sigma=~1)
par0 <- c(0.1,1,3)
response1<- SDE$new(formulas = formulas,data = dataAE[dataAE$ID="Asgeir',"],type = "CIR",
                    response = c("DistanceShore"),par0 = par0)



#fit_model
response1$fit(metho="Nelder-Mead")
estimates1=as.list(response1$tmb_rep(),what="Est")
std1=as.list(response1$tmb_rep(),what="Std")


#plot parameters
xmin=list("ExpShore"=1/1,"AngleNormal"=-pi+pi/20)
xmax=list("ExpShore"=1/0.05,"AngleNormal"=pi-pi/20)
link=list("ExpShore"=(\(x) 1/x))
xlabel=list("ExpShore"="Distance to shore")


res=response1$get_all_plots(baseline=NULL,model_name="response1",xmin=xmin,
                            xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)
