
library(sf)
library(ggplot2)

source("/home/delporta/Documents/Recherche/Codes/smoothSDE/R/utility.R")
# FJORDS POLYGONS

# Set the path to the directory containing coastline data
data_path <- file.path("/home","delporta", "Documents", "Recherche","Données", "greenland")

#get the coastline geometry from the geojson file
coastline<-st_read(file.path(data_path,"filtered_coastline_utm.geojson"))
land<-st_read(file.path(data_path,"filtered_land_utm.shp"))



#RECTANGULAR DOMAIN

#rectangular domain 
large_rectangle_coords <- matrix(c(0, 0, 0, 10000, 10000, 10000, 10000, 0, 0, 0), ncol = 2, byrow = TRUE)
small_rectangle_coords <- matrix(c(1000,1000, 1000, 9000, 9000, 9000, 9000, 1000, 1000, 1000), ncol = 2, byrow = TRUE)

# Create an sf object representing the rectangles
large_rectangle <- st_sfc(st_polygon(list(large_rectangle_coords)))
large_rectangle<-st_sf(geometry = large_rectangle)

small_rectangle <- st_sfc(st_polygon(list(small_rectangle_coords)))
small_rectangle<-st_sf(geometry = small_rectangle)

#change its crs
rect<- st_sym_difference(large_rectangle, small_rectangle)




# INITIAL POSITIONS

n_samples=6

#generate random initial points in the rectangle
z0_rect=matrix(rep(NA,n_samples*2),ncol=2)
colnames(z0_rect)=c("x1","x2")

for (i in 1:n_samples) {
  #choose location uniformly in the domain
  x=c(runif(1,min=1.5,max=8.5),runif(1,min=1.5,max=8.5))
  z0_rect[i,]=x
}

z0=matrix(as.numeric(c(dataAE[dataAE$ID=="Asgeir",][1,c("x","y")],dataAE[dataAE$ID=="Frederik",][1,c("x","y")],
                       dataAE[dataAE$ID=="Helge18",][1,c("x","y")],dataAE[dataAE$ID=="Kyrri",][1,c("x","y")],
                       dataAE[dataAE$ID=="Nemo",][1,c("x","y")],dataAE[dataAE$ID=="Siggi",][1,c("x","y")])),
          ncol=2,byrow=TRUE)



#FUNCTIONS TO COMPUTE COVARIATES ALONG THE WAY

fDshore=function(z,v,p) {
  
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  
  return (Dshore)
}

fExpShore=function(z,v,p) {
  
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  if (Dshore>0.5){
    return (0)
  }
  else if (Dshore<0.01) {
    return (100)
  }
  else {
    return (1/Dshore)
  }
}

fangle=function(z,v,p) {
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #return angle
  return (signed_angle(normal,v))
}


atw=list("ExpShore"=fExpShore,"AngleNormal"=fangle)


# TIME STEPS
Tmax=24
times=seq(0,Tmax,by=1/60)
n=length(times)
data_reg=data.frame("AngleNormal"=rep(0,n*length(unique(dataAE$ID))),
                    "ExpShore"=rep(1,n*length(unique(dataAE$ID))),times=seq(0,Tmax,by=1/60),ID=rep(unique(dataAE$ID),each=n))




# SIMULATE WITH CONSTRAINTS
sample1=response1$simulate(z0=z0,data=data_reg,atw=atw,land=land)
sample2=response2$simulate(z0=z0,data=data_reg,atw=atw,land=land)
sample1_rect=response1$simulate(z0=z0_rect,data=data_reg,atw=atw,land=rect)
sample2_rect=response2$simulate(z0=z0_rect,data=data_reg,atw=atw,land=rect)


#plots

p1=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample1,aes(x*1000,y*1000,col=ID),size=0.6)

p2=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample2,aes(x*1000,y*1000,col=ID),size=0.6)

p1_rect=ggplot()+geom_sf(data=rect,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample1_rect,aes(x*1000,y*1000,col=ID),size=0.6)

p2_rect=ggplot()+geom_sf(data=rect,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample2_rect,aes(x*1000,y*1000,col=ID),size=0.6)

saveWidget(ggplotly(p1), file ="response1/sample_response1.html")
saveWidget(ggplotly(p1_rect), file ="response1/sample_response1_rect.html")
saveWidget(ggplotly(p2), file ="response2/sample_response2.html")
saveWidget(ggplotly(p2_rect), file ="response2/sample_response2_rect.html")