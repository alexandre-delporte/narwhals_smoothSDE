

# FJORDS POLYGONS

# Set the path to the directory containing coastline data
data_path <- file.path("/home","delporta", "Documents", "Recherche","Données", "greenland")

#get the coastline geometry from the geojson file
coastline<-st_read(file.path(data_path,"filtered_coastline_utm.geojson"))
land<-st_read(file.path(data_path,"filtered_land_utm.shp"))



#RECTANGULAR DOMAIN
large_rectangle_coords <- matrix(c(0, 0, 0, 50000, 50000, 50000, 50000, 0, 0, 0), ncol = 2, byrow = TRUE)
small_rectangle_coords <- matrix(c(1000,1000, 1000, 49000, 49000, 49000, 49000, 1000, 1000, 1000), ncol = 2, byrow = TRUE)

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
  x=c(runif(1,min=1.5,max=48.5),runif(1,min=1.5,max=48.5))
  z0_rect[i,]=x
}

z0=matrix(as.numeric(c(dataAE[dataAE$ID=="Asgeir",][1,c("x","y")],dataAE[dataAE$ID=="Frederik",][1,c("x","y")],
                       dataAE[dataAE$ID=="Helge18",][1,c("x","y")],dataAE[dataAE$ID=="Kyrri",][1,c("x","y")],
                       dataAE[dataAE$ID=="Nemo",][1,c("x","y")],dataAE[dataAE$ID=="Siggi",][1,c("x","y")])),
          ncol=2,byrow=TRUE)



#FUNCTIONS TO COMPUTE COVARIATES ALONG THE WAY

fshore=function(z,v,p) {
  
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #distance to shore
  Dshore=sqrt(normal[1]^2+normal[2]^2)
  
  return (Dshore)
}

fangle=function(z,v,p) {
  
  #normal vector
  normal=c(z[1]-p[1],z[2]-p[2])
  
  #return angle
  return (signed_angle(normal,v))
}

atw=list("DistanceShore"=fshore,"AngleNormal"=fangle)


# TIME STEPS
Tmax=24
times=seq(0,Tmax,by=1/60)
n=length(times)
data_reg=data.frame("AngleNormal"=rep(0,6*n),"DistanceShore"=rep(1,6*n),times=seq(0,Tmax,by=1/60),ID=rep(unique(baseline0$data()$ID),each=n))



# SIMULATE WITHOUT CONSTRAINTS
sample1=response1$simulate(z0=z0,data=response1$data())
sample2=response2$simulate(z0=z0,data=response2$data())
sample3=response3$simulate(z0=z0,data=response3$data())
sample4=response4$simulate(z0=z0,data=response4$data())
sample5=response5$simulate(z0=z0,data=response5$data())
sample6=response6$simulate(z0=z0,data=response6$data())
sample7=response7$simulate(z0=z0,data=response7$data())
sample8=response8$simulate(z0=z0,data=response8$data())



# SIMULATE WITH CONSTRAINTS
sample9=response9$simulate(z0=z0,data=response9$data(),atw=atw,land=land)
sample10=response10$simulate(z0=z0,data=response10$data(),atw=atw,land=land)
sample9_rect=baseline9$simulate(z0=z0_rect,data=data_reg,atw=atw,land=rect)
sample10_rect=baseline10$simulate(z0=z0_rect,data=data_reg,atw=atw,land=rect)


#plots
p0=ggplot()+geom_point(data=dataAE,aes(x*1000,y*1000,col=ID),size=0.6)

p1=ggplot()+geom_point(data=sample1,aes(x*1000,y*1000,col=ID),size=0.6)

p2=ggplot()+geom_point(data=sample2,aes(x*1000,y*1000,col=ID),size=0.6)

p3=ggplot()+geom_point(data=sample3,aes(x*1000,y*1000,col=ID),size=0.6)

p4=ggplot()+geom_point(data=sample4,aes(x*1000,y*1000,col=ID),size=0.6)

p5=ggplot()+geom_point(data=sample5,aes(x*1000,y*1000,col=ID),size=0.6)

p6=ggplot()+geom_point(data=sample6,aes(x*1000,y*1000,col=ID),size=0.6)

p7=ggplot()+geom_point(data=sample7,aes(x*1000,y*1000,col=ID),size=0.6)

p8=ggplot()+geom_point(data=sample8,aes(x*1000,y*1000,col=ID),size=0.6)

p9=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample9,aes(x*1000,y*1000,col=ID),size=0.6)

p10=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample10,aes(x*1000,y*1000,col=ID),size=0.6)

p9_rect=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample9_rect,aes(x*1000,y*1000,col=ID),size=0.6)

p10_rect=ggplot()+geom_sf(data=land$geometry,fill=NA)+coord_sf(datum=st_crs(32626))+
  geom_point(data=sample10_rect,aes(x*1000,y*1000,col=ID),size=0.6)

saveWidget(ggplotly(p1), file ="response1/sample_response1.html")
saveWidget(ggplotly(p2), file ="response2/sample_response2.html")
saveWidget(ggplotly(p3), file ="response3/sample_response3.html")
saveWidget(ggplotly(p4), file ="response4/sample_response4.html")
saveWidget(ggplotly(p5), file ="response5/sample_response5.html")
saveWidget(ggplotly(p6), file ="response6/sample_response6.html")
saveWidget(ggplotly(p7), file ="response7/sample_response7.html")
saveWidget(ggplotly(p8), file ="response8/sample_response8.html")
saveWidget(ggplotly(p9), file ="response9/sample_response9.html")
saveWidget(ggplotly(p9_rect), file ="response9/sample_response9_rect.html")
saveWidget(ggplotly(p10), file ="response10/sample_response10.html")
saveWidget(ggplotly(p10_rect), file ="response10/sample_response10_rect.html")