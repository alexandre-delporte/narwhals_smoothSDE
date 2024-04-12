

#SEED FOR REPRODUCTIBILITY
set.seed(42)

#REQUIRED LIBRARIES
library(smoothSDE)
library(ggplot2)
#source(file.path("/home", "delporta","Documents","Recherche","Codes","smoothSDE","R","utility.R"))



################################### OUTPUT FILE TO KEEP RECORD OF THE SCRIPT EXECUTION ####################
filename=paste("fit_response_output.txt",sep="")

#Delete file if it exists
if (file.exists(filename)) {
  file.remove(filename)
}

file.create(filename)

########################################    GET NARWHAL DATA   #############################################
par_dir=dirname(dirname(dirname(getwd())))
narwhals_data_path<- file.path(par_dir,"Data","Narwhals")       

# DATA BEFORE EXPOSURE

dataBE1=read.csv(file.path(narwhals_data_path,"DataBE1.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 1 day after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE1[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)

dataBE2=read.csv(file.path(narwhals_data_path,"DataBE2.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories before exposure and 12h after tagging... ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataBE2[,1]),"positions measured before exposure \n"),file=filename,sep="\n",append=TRUE)


# DATA AFTER EXPOSURE

dataAE=read.csv(file.path(narwhals_data_path,"DataAE.csv"), header = TRUE,dec = ".")

cat("Extracting trajectories after exposure.. ",file=filename,sep="\n",append=TRUE)

cat(paste(length(dataAE[,1]),"positions measured after exposure \n"),file=filename,sep="\n",append=TRUE)




############################################# SET MEASUREMENT ERROR ###############################################

sigma_obs=0.03
n_obs=length(dataAE$time)
H=array(rep(sigma_obs^2*diag(2),n_obs),dim=c(2,2,n_obs))

####################################  RESPONSE MODELS WITH ExpShip IN NU ###########################################


#initialize parameters with estimated baseline values
par0=sapply(1:4, function(i) {baseline2$invlink()[[i]](estimates2$coeff_fe[i])})

#Find best degree of freedom for cubic splines

list_k=c(3,4,7,10,15,20,30)
AIC_df1=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {
  
  formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re"),tau = ~1,
                   nu=~s(ExpShip,k=list_k[i],bs="cs")+s(ID,bs="re"))

  try_response1 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                       par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
  try_response1$fit()
  
  AIC_df1[i,"AIC_marginal"]=try_response1$AIC_marginal()
  AIC_df1[i,"AIC_conditional"]=try_response1$AIC_conditional()
}

plotAIC1=ggplot(data=AIC_df1,aes(x=k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(k, AIC_conditional,color="AIC_conditional"))+geom_point(aes(k, AIC_conditional),shape=21)+
  xlab("Cubic splines degree of freedom")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)
  
ggsave(path="response1","plotAIC1.png",plot=plotAIC1)


#MODEL WITHOUT OFFSET


formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re"),tau = ~1,
                 nu=~s(ExpShip,k=6,bs="cs")+s(ID,bs="re"))
response1 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response1$fit()


xmin=list("ExpShip"=1/80)
xmax=list("ExpShip"=1/3)
link=list("ExpShip"=(\(x) 1/x))
xlabel=list("ExpShip"="Distance to ship")

plot_resp1=response1$get_all_plots(baseline=baseline2,model_name="response1",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


## MODEL WITH OFFSET AS IN BASELINE 2

map=list("coeff_fe"=factor(rep(NA,4)))

response1offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",
                           response = c("x","y"),par0 = par0,map=map,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response1offset$fit()

plot_resp1offset=response1offset$get_all_plots(baseline=baseline2,model_name="response1offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


##################################  RESPONSE MODELS WITH ExpShip IN TAU  #################################################


#Find best degree of freedom for cubic splines

AIC_df2=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {

  formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau =~s(ExpShip,k=list_k[i],bs="cs"),
                   nu=~s(ID,bs="re"))
  
  try_response2 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
  try_response2$fit()
  
  AIC_df2[i,"AIC_marginal"]=try_response2$AIC_marginal()
  AIC_df2[i,"AIC_conditional"]=try_response2$AIC_conditional()
}

plotAIC2=ggplot(data=AIC_df2,aes(x=k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(k, AIC_conditional,color="AIC_conditional"))+geom_point(aes(k, AIC_conditional),shape=21)+
  xlab("Cubic splines degree of freedom")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)

ggsave(path="response2","plotAIC2.png",plot=plotAIC2)


#set k to best degree of freedom
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau =~s(ExpShip,k=7,bs="cs"),
                 nu=~s(ID,bs="re"))

## FIT MODEL WITHOUT OFFSET

response2 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response2$fit()

#get parameters plots
plot_resp2=response2$get_all_plots(baseline=baseline2,model_name="response2",xmin=xmin,xmax=xmax,
              link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)

## FIT MODEL WITH OFFSET AS IN BASELINE 1
response2offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),map=map)
response2offset$fit()

#get parameters plot
plot_resp2offset=response2offset$get_all_plots(baseline=baseline2,model_name="response2offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


###################################    RESPONSE MODELS WITH ExpShip IN TAU AND NU ##############################################



#Find best degree of freedom for cubic splines

AIC_df3=data.frame("k"=list_k,"AIC_marginal"=rep(NA,length(list_k)),"AIC_conditional"=rep(NA,length(list_k)))

for (i in 1:length(list_k)) {
  
  
  formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~s(ExpShip,k=list_k[i],bs="cs"),
                   nu=~s(ExpShip,k=list_k[i],bs="cs")+s(ID,bs="re"))
  
  try_response3 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
  try_response3$fit()
  
  AIC_df3[i,"AIC_marginal"]=try_response3$AIC_marginal()
  AIC_df3[i,"AIC_conditional"]=try_response3$AIC_conditional()
}

plotAIC3=ggplot(data=AIC_df3,aes(x=2*k))+geom_line(aes(y=AIC_marginal,color="AIC_marginal"))+geom_point(aes(y=AIC_marginal),shape=21)+
  geom_line(aes(y=AIC_conditional,color="AIC_conditional"))+geom_point(aes(y=AIC_conditional),shape=21)+
  xlab("Number of coefficients to estimate (=2q)")+ylab("AIC")+
  scale_color_manual(values = c("AIC_marginal" = "red", "AIC_conditional" = "blue"),
                     labels = c("Marginal AIC", "Conditional AIC"),name=NULL)

ggsave(path="response3","plotAIC3.png",plot=plotAIC3,width=10,height=5)


#Set k to best degree of freedom
formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re"),tau = ~s(ExpShip,k=10,bs="cs"),
                 nu=~s(ExpShip,k=10,bs="cs"))



## FIT MODEL WITHOUT OFFSET
response3 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response3$fit()

plot_resp3=response3$get_all_plots(baseline=baseline2,model_name="response3",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)

#inspect specific values of tau
tau_data=plot_resp3$data_tau_ExpShip
q1=quantile(dataAE$Dist_Ship,0.05,na.rm=TRUE)
q2=quantile(dataAE$Dist_Ship,0.25,na.rm=TRUE)
q3=quantile(dataAE$Dist_Ship,0.75,na.rm=TRUE)
very_close_tau=mean(tau_data[1/tau_data$ExpShip<q1,"par"])
close_tau=mean(tau_data[1/tau_data$ExpShip>q1 & 1/tau_data$ExpShip<q2,"par"])
medium_tau=mean(tau_data[1/tau_data$ExpShip>q2 & 1/tau_data$ExpShip<q3,"par"])
far_tau=mean(tau_data[1/tau_data$ExpShip>q3,"par"])

#inspect specific values of nu
nu_data=plot_resp3$data_nu_ExpShip
very_close_nu=mean(nu_data[1/nu_data$ExpShip<q1,"par"])
close_nu=mean(nu_data[1/nu_data$ExpShip>q1 & 1/nu_data$ExpShip<q2,"par"])
medium_nu=mean(nu_data[1/nu_data$ExpShip>q2 & 1/nu_data$ExpShip<q3,"par"])
far_nu=mean(nu_data[1/nu_data$ExpShip>q3,"par"])

## FIT MODEL WITH OFFSET AS IN BASELINE 2
response3offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response3offset$fit()

plot_resp3offset=response3offset$get_all_plots(baseline=baseline2,model_name="response3offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


###############################   RESPONSE MODELS WITH XT AND XI COVARIATE IN NU #############################################

formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~1,
                 nu=~s(ExpShipT,k=7,bs="cs")+
                   s(ExpShipI,k=7,bs="cs")+s(ID,bs="re"))

## MODEL WITHOUT OFFSET
response4 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response4$fit()


xmin=list("ExpShipT"=1/80,"ExpShipI"=1/80)
xmax=list("ExpShipT"=1/3,"ExpShipI"=1/3)
link=list("ExpShipT"=(\(x) 1/x),"ExpShipI"=(\(x) 1/x))
xlabel=list("ExpShipT"="Distance to ship","ExpShipI"="Distance to ship")

plot_resp4=response4$get_all_plots(baseline=baseline2,model_name="response4",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


## MODEL WITH OFFSET AS IN BASELINE 1

response4offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response4offset$fit()

plot_resp4offset=response4offset$get_all_plots(baseline=baseline1,model_name="response4offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


###################################  RESPONSE MODEL WITH ExpShipT AND ExpShipI COVARIATES IN TAU  ####################################


formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~
                   s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs"),
                 nu=~s(ID,bs="re"))

## MODEL WITHOUT OFFSET
response5 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response5$fit()


plot_resp5=response5$get_all_plots(baseline=baseline2,model_name="response5",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


## MODEL WITH OFFSET AS IN BASELINE 1

response5offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)),map=map)
response5offset$fit()


plot_resp5offset=response5offset$get_all_plots(baseline=baseline2,model_name="response5offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)



#############################   RESPONSE  MODEL WITH ExpShipT AND ExpShipI COVARIATES IN NU AND TAU    ####################################


formulas <- list(mu1 = ~s(ID,bs="re") ,mu2 =~s(ID,bs="re") ,tau = ~s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs")
                 ,nu=~s(ExpShipT,k=7,bs="cs")+s(ExpShipI,k=7,bs="cs")+s(ID,bs="re"))

## MODEL WITHOUT OFFSET

response6 <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                     par0 = par0,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response6$fit()


plot_resp6=response6$get_all_plots(baseline=baseline2,model_name="response6",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


## MODEL WITH OFFSET AS IN BASELINE 1

response6offset <- SDE$new(formulas = formulas,data = dataAE,type = "CTCRW",response = c("x","y"),
                           par0 = par0,map=map,other_data=list("log_sigma_obs0"=log(sigma_obs)))
response6offset$fit()


plot_resp6offset=response6offset$get_all_plots(baseline=baseline2,model_name="response6offset",
              xmin=xmin,xmax=xmax,link=link,xlabel=xlabel,show_CI="pointwise",save=TRUE)


#############################      RESPONSE AIC VALUES    #####################################################

AIC_resp1=response1$AIC_marginal()
AIC_resp2=response2$AIC_marginal()
AIC_resp3=response3$AIC_marginal()
AIC_resp4=response4$AIC_marginal()
AIC_resp5=response5$AIC_marginal()
AIC_resp6=response6$AIC_marginal()



AIC_resp1offset=response1offset$AIC_marginal()
AIC_resp2offset=response2offset$AIC_marginal()
AIC_resp3offset=response3offset$AIC_marginal()
AIC_resp4offset=response4offset$AIC_marginal()
AIC_resp5offset=response5offset$AIC_marginal()
AIC_resp6offset=response6offset$AIC_marginal()






cat("--------------------- RESPONSE MODELS BASED ON BASELINE 2 ---------------------: \n",
    "WITHOUT OFFSET : \n",response1$message(),AIC_resp1,
    response2$message(),AIC_resp2,response3$message(),AIC_resp3,
    response4$message(),AIC_resp4,response5$message(),AIC_resp5,
    response6$message,AIC_resp6,file=filename,sep="\n",append=TRUE)

cat("\n",file=filename,append=TRUE)


cat("WITH OFFSET : \n",response1offset$message(),AIC_resp1offset,
    response2offset$message(),AIC_resp2offset,response3offset$message(),AIC_resp3offset,
    response4offset$message(),AIC_resp4offset,response5offset$message(),AIC_resp5offset,
    response6offset$message,AIC_resp6offset,file=filename,sep="\n",append=TRUE)


