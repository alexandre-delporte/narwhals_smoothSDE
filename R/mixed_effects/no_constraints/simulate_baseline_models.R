

library(dplyr)
library(smoothSDE)
library(plotly)
library(htmlwidgets)


####################################    GET NARWHAL DATA   ###########################################

# Set the path to the directory containing the data
par_dir=dirname(dirname(dirname(getwd()))) #parent directory 
narwhal_data_path <- file.path(par_dir,"Data", "Narwhals")

# DATA BEFORE EXPOSURE
dataBE1=read.csv(file.path(narwhal_data_path,"DataBE1.csv"), header = TRUE,dec = ".")


# INITIAL POSITIONS
z0=matrix(as.numeric(c(dataBE1[dataBE1$ID=="Asgeir",][1,c("x","y")],dataBE1[dataBE1$ID=="Frederik",][1,c("x","y")],
          dataBE1[dataBE1$ID=="Helge18",][1,c("x","y")],dataBE1[dataBE1$ID=="Kyrri",][1,c("x","y")],
          dataBE1[dataBE1$ID=="Nemo",][1,c("x","y")],dataBE1[dataBE1$ID=="Siggi",][1,c("x","y")])),
          ncol=2,byrow=TRUE)


# TIME STEPS
Tmax=5*24
times=seq(0,Tmax,by=3/60)
n=length(times)
data_reg=data.frame(time=times,ID=rep(unique(baseline0$data()$ID),each=n))


# SIMULATE 
sample0=baseline0$simulate(z0=z0,data=data_reg)
sample1=baseline1$simulate(z0=z0,data=data_reg)
sample2=baseline2$simulate(z0=z0,data=data_reg)
sample3=baseline3$simulate(z0=z0,data=data_reg)
sample4=baseline4$simulate(z0=z0,data=data_reg)

#plots
pBE=ggplot()+geom_point(data=dataBE1,aes(x*1000,y*1000,col=ID),size=0.6)

p0=ggplot()+geom_point(data=sample0,aes(x*1000,y*1000,col=ID),size=0.6)+
  geom_point(data = sample0 %>% filter(!duplicated(ID)),
             aes(x = x * 1000, y = y * 1000), shape = 3, size = 3, col = "black")

p1=ggplot()+geom_point(data=sample1,aes(x*1000,y*1000,col=ID),size=0.6)+
  geom_point(data = sample1 %>% filter(!duplicated(ID)),
             aes(x = x * 1000, y = y * 1000), shape = 3, size = 3, col = "black")

p2=ggplot()+geom_point(data=sample2,aes(x*1000,y*1000,col=ID),size=0.6)+
  geom_point(data = sample2 %>% filter(!duplicated(ID)),
             aes(x = x * 1000, y = y * 1000), shape = 3, size = 3, col = "black")

p3=ggplot()+geom_point(data=sample3,aes(x*1000,y*1000,col=ID),size=0.6)+
  geom_point(data = sample3 %>% filter(!duplicated(ID)),
             aes(x = x * 1000, y = y * 1000), shape = 3, size = 3, col = "black")
  
p4=ggplot()+geom_point(data=sample4,aes(x*1000,y*1000,col=ID),size=0.6)+
  geom_point(data = sample4 %>% filter(!duplicated(ID)),
             aes(x = x * 1000, y = y * 1000), shape = 3, size = 3, col = "black")

saveWidget(ggplotly(p1), file ="baseline1/sample_baseline1.html")
saveWidget(ggplotly(p2), file ="baseline2/sample_baseline2.html")
saveWidget(ggplotly(p3), file ="baseline3/sample_baseline3.html")
saveWidget(ggplotly(p4), file ="baseline4/sample_baseline4.html")


