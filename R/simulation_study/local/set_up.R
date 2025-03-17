

library(here)
library(sf)
library(smoothSDE)
library(doParallel)
library(ggplot2)
library(dplyr)
library(mgcv)
library(plotly)

#get functions to simulate trajectories
source(file.path(here("R","simulation_study","CVM_functions.R")))  


## Define hyperparameters 

N_ID=6

TMAX=24
DELTA=1/60

SP_DF=c(5,5)

TAU_0=1.5
NU_0=4

SIGMA_TAU=0.2
SIGMA_NU=0.1

DMIN=0.5
DMAX=1

PAR0_RACVM=c(0,0,1,1,0)

SIGMA_OBS_LOW=0.01
SIGMA_OBS_HIGH=0.04

BY_LF=5
BY_HF=1

D_LOW=0.05
D_UP=3

A=2
B=3
D0=0.5
D1=1
SIGMA_THETA=pi/6
SIGMA_D=0.3



seed=2
set.seed(seed)
# Set the path to the directory containing the data
greenland_data_path <- here("Data","preprocessed_data","greenland")

#get the land geometry from shapefile
border<-st_read(file.path(greenland_data_path,"updated_scoresby_sound_utm.shp"))
border <- st_transform(border, crs = "+init=EPSG:32626 +units=km")

v0=c(0,0)

#generate random initial points
x0=matrix(rep(NA,N_ID*2),ncol=2)
colnames(x0)=c("x1","x2")

i=1
while (i<=N_ID) {
  #choose location uniformly in the map
  x=c(runif(1,min=430,max=500),runif(1,min=7760,max=7900))
  if (!(is_in_border(st_point(x),border))) {
    p=nearest_boundary_points(matrix(x,ncol=2),border)
    Dshore=(x[1]-p[1])^2+(x[2]-p[2])^2
    #keep it as initial position if it is at least 50 metres away from the shore
    if (Dshore>DMIN && Dshore<DMAX) {
      x0[i,]=x
      i=i+1
    }
  }
}

#high frequency time points
times=seq(0,TMAX,by=DELTA)

n_obs=length(times)-1

# Definition of parameters tau and nu ----------------
tau_re=rnorm(N_ID,mean=0,sd=SIGMA_TAU)
nu_re=rnorm(N_ID,mean=0,sd=SIGMA_NU)
true_log_tau=tau_re+log(TAU_0)
true_log_nu=nu_re+log(NU_0)


# Defintion of smooth parameter omega ----------


fomega_cubic=function(cov_data,boundary_distance_name="BoundaryDistance",
                      boundary_angle_name="BoundaryAngle",
                      a,b,D0,D1,sigma_D,sigma_theta){
  Dshore=cov_data[[boundary_distance_name]]
  theta=cov_data[[boundary_angle_name]]
  
  a*theta*(theta-pi/2)*(theta+pi/2)*exp(-Dshore/D0)/Dshore+
    b*(exp(-1/2*(((theta+pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2))-
         exp(-1/2*(((theta-pi/2/sqrt(3))/sigma_theta)^2+((Dshore-D1)/sigma_D)^2)))
}


fomega=function(cov_data,boundary_distance_name="BoundaryDistance",
                boundary_angle_name="BoundaryAngle") {
  fomega_cubic(cov_data,boundary_distance_name,
               boundary_angle_name,
               A,B,D0,D1,SIGMA_D,SIGMA_THETA)
}


# Generate samples ---------------

cores=detectCores()
cl <- makeCluster(6) #not to overload your computer
registerDoParallel(cl)



data=foreach (i=1:N_ID,.combine='rbind',.packages=c("progress","MASS","sf","mgcv","smoothSDE")) %dopar% {
  
  set.seed((seed-1)*N_ID+i)
  
  #constant nu
  fnu_constant=function(cov_data) {
    return (exp(true_log_nu[i]))
  }
  
  ftau_constant=function(cov_data) {
    return (exp(true_log_tau[i]))
  }
  
  
  res=sim_CRCVM(ftau=ftau_constant,fomega=fomega,fnu=fnu_constant,
                log_sigma_obs=NULL,v0=v0,x0=x0[i,],times=times,
                land=border,verbose=FALSE,save_omega = TRUE)
  
  data_sim=res$sim
  data_sim$ID=factor(rep(i,length(data_sim$y1)))
  data_shore=res$boundary
  
  data_sim=cbind(data_sim,data_shore)
  data_sim
}

#stop cluster
stopCluster(cl)




# Points that reached land ---------------
count=0
remove_ID=c()
for (id in unique(data$ID)) {
  id_data=data[data$ID==id,]
  if (nrow(id_data) < n_obs) {
    count=count+1
    remove_ID=c(remove_ID,id)
    cat("ID",id,"reached land","\n",sep=" ")
  }
}

cat(count/N_ID*100,"percent of the trajectories reached land","\n")


data_sub=data %>%  filter(!(ID %in% remove_ID))


# add noise
low_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_LOW^2,2)))
high_noise=rmvn(n_obs,rep(0,2),diag(rep(SIGMA_OBS_HIGH^2,2)))
observed_data_le=data_sub[,c("y1","y2","time","ID","true_BoundaryDistance","true_BoundaryAngle")]
observed_data_le[,c("y1","y2")]=data_sub[,c("y1","y2")]+low_noise
observed_data_he=data_sub[,c("y1","y2","time","ID","true_BoundaryDistance","true_BoundaryAngle")]
observed_data_he[,c("y1","y2")]=data_sub[,c("y1","y2")]+high_noise

#subsample
data_lf_he=observed_data_he[seq(1,length(data_sub$time),by=BY_LF),]
data_lf_le=observed_data_le[seq(1,length(data_sub$time),by=BY_LF),]
data_hf_le=observed_data_le[seq(1,length(data_sub$time),by=BY_HF),]
data_hf_he=observed_data_he[seq(1,length(data_sub$time),by=BY_HF),]


# Save plot of the trajectories ------------

plot=ggplot()+geom_sf(data=border$geometry,fill="grey")+
  coord_sf(datum=st_crs("+init=EPSG:32626 +units=km"))+
  geom_point(data=data_sub,mapping=aes(y1,y2,col=ID),size=0.1,alpha=0.2,shape=2)+
  geom_point(data=data_lf_he,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_path(data=data_lf_he,mapping=aes(y1,y2,col=ID),size=0.1)+
  geom_point(data = data%>% filter(!duplicated(ID)),
             aes(x = y1, y = y2), shape = 3, size = 4, col = "red")+
  xlab("x") + ylab("y")

ggsave(filename="plot_trajectories.pdf",plot=plot,width=10,height=5) 

ggplotly(plot)


data_list=list("data_lf_he"=data_lf_he,"data_hf_le"=data_hf_le,"data_lf_le"=data_lf_le,
               "data_hf_he"=data_hf_he)

data_list=lapply(data_list,function(data) {
  covs=get_BoundaryMetrics(data,response=c("y1","y2"),border=border,n_cores=6)
  colnames(covs)=c("observed_BoundaryDistance","observed_BoundaryAngle")
  data=cbind(data,covs)
})


plot_Dshore=ggplot(data=data_list[["data_lf_he"]])+
  geom_point(aes(x=time,y=observed_BoundaryDistance,col=ID),size=1,alpha=0.5)+
  geom_line(aes(x=time,y=observed_BoundaryDistance,col=ID),size=1,alpha=0.2)+
  facet_wrap(~ID)+
  theme_minimal()

hist_Dshore=ggplot(data=data_list[["data_lf_he"]])+
  geom_histogram(aes(x=observed_BoundaryDistance),bins=30)+
  facet_wrap(~ID)+theme_minimal()

plot_Dshore
hist_Dshore
