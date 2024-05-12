
#SCRIPT CONTAINING FUNCTIONS TO SIMULATE AND COMPUTE LIKELIHOOD FOR DIFFERENT CVM MODELS 

library(Matrix)
library(MASS)
library(progress)




#' @description Compute likelihood for a 2D CVM with constant parameters and Kalman filter
#' 
#' @param obs data frame with at least columns matching with the names in the response argument 
#' and one "time" column
#' @param response vector of strings corresponding to the names of the response
#' @param log_beta log of autocorrelation time scale parameter
#' @param log_sigma log of of volatility parameter
#' @param mu1 first component of mean velocity: scalar with size length(times)-1nent of the drift parameter,
#'  scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of mean velocity: scalar with size length(times)-1nent of the drift parameter, 
#' scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param log_sigma_obs if not NULL, gaussian noise with standard deviation exp(log_sigma_obs) is added to the sampled locations 
#' @return scalar corresponding to the likelihood of the data

nllkCVM=function(obs,response,log_beta,log_sigma,mu1,mu2,v0,x0,log_sigma_obs) {
  
  
  n=length(obs[,response[1]])
  
  beta=exp(log_beta)
  sigma=exp(log_sigma)
  mu=matrix(c(mu1,mu2),ncol=1,nrow=2,byrow=TRUE)
  
  times=obs$time
  
  
  #initialization
  alphax=c(x0[1],v0[1]-mu1)
  alphay=c(x0[2],v0[2]-mu2)
  alpha=c(alphax,alphay)
  Rx=matrix(rep(0,4),ncol=2)
  Ry=matrix(rep(0,4),ncol=2)
  R=bdiag(Rx,Ry)
  #log likelihood
  llk=-n*log(2*pi)
  
  for (i in (2:n)) {
    #time step
    delta=times[i]-times[i-1]
    
    #link matrix for successive states
    Ti=matrix(c(1,1/beta*(1-exp(-beta*delta)),0,exp(-beta*delta)),nrow=2,byrow=TRUE)
    Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
    
    #observation - states link matrix 
    Z_tilde=matrix(c(1,0,0,0,0,0,1,0),ncol=4,byrow=TRUE)
    
    #covariance matrix eta
    var_xi=sigma^2/beta^2*(delta+2/beta*(exp(-beta*delta)-1)+1/2/beta*(1-exp(-2*beta*delta)))
    var_zeta=sigma^2/2/beta*(1-exp(-2*beta*delta))
    cov=sigma^2/2/beta^2*(1+exp(-2*beta*delta)-2*exp(-beta*delta))
    Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
    Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
    
    #drift part
    Bi=matrix(c(delta,0),nrow=2,ncol=1,byrow=TRUE)
    Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
    
    #measurement error matrix
    if (!(is.null(log_sigma_obs))) {
      Hi=diag(exp(log_sigma_obs),2)
    }
    else {
      Hi=matrix(c(0,0,0,0),ncol=2)
    }
    
    #observation 
    obsi=matrix(c(obs[i,response[1]],obs[i,response[2]]),nrow=2,ncol=1,byrow=TRUE)
    
    
    #prediction
    alpha_=Ti_tilde%*%alpha+Bi_tilde%*%mu
    R_=Ti_tilde%*%R%*%t(Ti_tilde)+Qi_tilde
    
    #correction
    Fi=Z_tilde%*%R_%*%t(Z_tilde)+Hi
    Fiinv=solve(Fi)
    
    K=R_%*%t(Z_tilde)%*%Fiinv
    ui=obsi-Z_tilde%*%alpha_
    alpha=alpha_+K%*%(ui)
    R=(diag(4)-K%*%Z_tilde)%*%R_
    #update loglikelihood
    llk=llk-1/2*(log(abs(det(Fi)))+t(ui)%*%Fiinv%*%ui)
  }
  
  return (as.numeric(-llk))
  
}

#' @description Simulate from a 2D CVM
#' 
#' @param mu1 first component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param times vector of sampling times
#' @param log_sigma_obs if not NULL, gaussian noise with standard deviation exp(log_sigma_obs) is added to the sampled locations 
#' @param keep_true_pos default FALSE. If TRUE, return observations (with errors) and true positions in the dataframa
#' @param verbose if TRUE, display computations along the way.
#' @return dataframe with dim (n,3) where n is the length of 'times' argument and columns "y1", "y2" are the sampled
#' observations (with or without noise) and "time" is the time of the samples if keep_true_pos is FALSE, *
#' or (n,5) and the two additional columns "x1", "x2" are the trues positions (without noise)

sim_CVM=function(mu1,mu2,beta,sigma,v0,x0,times,log_sigma_obs=NULL,verbose=FALSE,keep_true_pos=FALSE) {
  
  
  n=length(times)
  
  betas=beta
  mu1s=mu1
  mu2s=mu2
  sigmas=sigma
  
  if (length(mu1)==1 & length(mu2)==1) {
    mu1s=rep(mu1,n)
    mu2s=rep(mu2,n)
  }
  if (length(beta)==1) {
    betas=rep(beta,n)
  }
  if (length(sigma)==1) {
    sigmas=rep(sigma,n)
  }
  
  par=data.frame(beta=betas,sigma=sigmas,mu1=mu1s,mu2=mu2s)
  
  mu=matrix(c(par$mu1[1],par$mu2[1]),ncol=1,nrow=2,byrow=TRUE)
 
  alpha=matrix(c(x0[1],v0[1]-mu[1],x0[2],v0[2]-mu[2]),nrow=4,ncol=1,byrow=TRUE)
  
  if (is.null(log_sigma_obs)){
    yobs=x0
  }
  else {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  
  data_sim=data.frame("y1"=yobs[1],"y2"=yobs[2],"x1"=alpha[1],"x2"=alpha[3])
  
  if (!(verbose)) {
    
    #progress bar
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  #loop over time steps
  for (i in 2:n) {
    
    if (!(verbose)) {
      pb$tick()
    }
    else {
      cat("velocity :",c(alpha[2]+mu1s[i-1],alpha[4]+mu2s[i-1]))
      cat("position :",c(alpha[1],alpha[3]),"\n")
    }
    
    #time step
    delta=times[i]-times[i-1]
    
    #parameters
    beta=par$beta[i-1]
    sigma=par$sigma[i-1]
    mu=matrix(c(par$mu1[i-1],par$mu2[i-1]),ncol=1,nrow=2,byrow=TRUE)
    
    #covariance matrix
    var_xi=sigma^2/beta^2*(delta+2/beta*(exp(-beta*delta)-1)+1/2/beta*(1-exp(-2*beta*delta)))
    var_zeta=sigma^2/2/beta*(1-exp(-2*beta*delta))
    cov=sigma^2/2/beta^2*(1+exp(-2*beta*delta)-2*exp(-beta*delta))
    Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
    Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
  
    #drift part
    Bi=matrix(c(delta,0),nrow=2,ncol=1,byrow=TRUE)
    Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
    
    #random part
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi_tilde)
    
    #hidden states link matrix
    Ti=matrix(c(1,1/beta*(1-exp(-beta*delta)),0,exp(-beta*delta)),nrow=2,byrow=TRUE)
    Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
    
    alpha=Ti_tilde%*%alpha+Bi_tilde%*%mu+eta
    
    if (is.null(log_sigma_obs)){
      yobs=c(alpha[1],alpha[3])
    }
    else {
      #measurement error
      epsilon=mvrnorm(2,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #state obs link matrix
      Z=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
      Z_tilde=bdiag(replicate(2,Z,simplify=FALSE))
      
      yobs=Z_tilde%*%alpha+epsilon
    }
   
    data_sim=rbind(data_sim,c(yobs[1],yobs[2],alpha[1],alpha[3]))
  }

  data_sim$time=times
  
  if (keep_true_pos) {
    return (data_sim)}
  else {
    return (data_sim[,c("y1","y2","time")])
  }
}



#' @description Compute nearest shore point 
#' 
#' @param point the location from which we want to compute the nearest shore point (as given by st_point)
#' @param coastline sf object (list of polygons) defining the land or the coastline
#' @return vector of coordinates of the nearest point on the shore

nearest_shore_point=function(point,coastline) {
  
  # Initialize variables to store the minimum distance and nearest point
  min_distance <- Inf
  nearest_point <- NULL
  
  
  # Iterate over each polygon in coastline
  for (j in 1:length(coastline$geometry)) {
    
    # Get the current polygon
    polygon <- coastline$geometry[[j]]
    
    
    distance <- st_distance(point, polygon)
    
    candidate=st_nearest_points(point,polygon)
    
    if (!st_is_empty(candidate)) {
    
      coordinates=st_coordinates(candidate)[2,c("X","Y")]
    
      if (distance < min_distance) {
        min_distance=distance
        nearest_point=coordinates
      }
    }
  }
  return (as.matrix(nearest_point))
}


#' @description Check whether a point is on land or not
#' 
#' @param point the location we want to tell if it is on land or not
#' @param land sf object (list of polygons) defining the land
#' @return boolean TRUE if it is on land, FALSE otherwise

is_in_land=function(point,land) {
  return (length(st_intersects(point,land$geometry)[[1]])>0)
}



#' @description Compute distances to shore 
#' 
#' @param df dataframe of GPS positions in UTM and km
#' @param response columns of the datafram containing the GPS positions
#' @param land sf object (list of polygons) defining the land (in UTM and metres)
#' @return vector of distances to shore

get_shore_distances=function(df,response,land) {
  
  positions=df[,response]
  n=length(positions[,1])
  Dshore=c()
  
  #progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent",
    total = n, clear = FALSE, width = 60)
  
  
  
  for (i in 1:n) {
    
    pb$tick()
    
    #position in km
    pos=as.numeric(positions[i,])
    
    #project
    proj=c(nearest_shore_point(st_point(1000*pos),land)/1000)
    
    #distance
    Dshore=c(Dshore,sqrt((pos[1]-proj[1])^2+(pos[2]-proj[2])^2))
  }
  return (Dshore)
}


#' @description Compute signed angle in [-pi,pi] that rotates first vector into second vector as in 
#' https://math.stackexchange.com/questions/529555/signed-angle-between-2-vectors
#' 
#' @param u vector or matrice of vectors of two columns
#' @param v vector matrice of bivariate vectors with 2 columns
#' @return angle in [-pi,pi]
#' 
signed_angle <- function(u, v) {
  
  u <- matrix(u, ncol = 2)
  v <- matrix(v, ncol = 2)
  if (nrow(u) != nrow(v)) stop("u and v must have the same number of rows")
  result <- as.numeric(atan2(v[,2], v[,1]) - atan2(u[,2], u[,1]))
  ind1 <- which(result > pi)
  ind2 <- which(result <= -pi)
  result[ind1] <- result[ind1] - 2*pi
  result[ind2] <- result[ind2] + 2*pi
  return(result) 
} 

#' @description Simulate from a 2D reflected CVM (projection when process gets on the land)
#' 
#' @param mu1 first component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param times vector of sampling times
#' @param log_sigma_obs if not NULL, gaussian noise with standard deviation exp(log_sigma_obs) is added to the sampled locations 
#' @param land sf object (list of polygons) defining the land
#' @param verbose if TRUE, display computations along the way.
#' @return dataframe with dim (n,3) where n is the length of 'times' argument and columns "y1", "y2" are the sampled
#' observations (with or without noise) and "time" is the time of the samples.

sim_reflectedCVM=function(mu1,mu2,beta,sigma,v0,x0,times,log_sigma_obs=NULL,land,coastline,verbose=TRUE) {
  
  #number of measurements
  n=length(times)
  
  betas=beta
  sigmas=sigma
  mu1s=mu1
  mu2s=mu2
 
  if (length(beta)==1) {
    betas=rep(beta,n-1)
  }
  if (length(sigma)==1) {
    sigmas=rep(sigma,n-1)
  }
  if (length(mu1)==1) {
    mu1s=rep(mu1s,n-1)
  }
  if (length(mu2)==1) {
    mu2s=rep(mu2s,n-1)
  }
  
  alpha=matrix(c(x0[1],v0[1]-mu1s[1],x0[2],v0[2]-mu2s[1]),nrow=4,ncol=1,byrow=TRUE)
  
  if (is.null(log_sigma_obs)) {
    yobs=c(alpha[1],alpha[3])
  }
  else {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
  
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  
  data_sim=data.frame("y1"=yobs[1],"y2"=yobs[2])
  
  
  if (!(verbose)) {
  #progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent",
    total = n, clear = FALSE, width = 60)
  }
  
  for (i in 2:n) {
    
    if (!(verbose)) {
      pb$tick()
    }
  
    #time step
    delta=times[i]-times[i-1]
    
    #previous position in metres
    prev_pos=1000*c(data_sim$x1[i-1],data_sim$x2[i-1])
    
    #parameters
    beta=betas[i-1]
    sigma=sigmas[i-1]
    mu1=mu1s[i-1]
    mu2=mu2s[i-1]
    mu=matrix(c(mu1,mu2),ncol=1,nrow=2,byrow=TRUE)
    
    #covariance matrix
    var_xi=sigma^2/beta^2*(delta+2/beta*(exp(-beta*delta)-1)+1/2/beta*(1-exp(-2*beta*delta)))
    var_zeta=sigma^2/2/beta*(1-exp(-2*beta*delta))
    cov=sigma^2/2/beta^2*(1+exp(-2*beta*delta)-2*exp(-beta*delta))
    Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
    Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
    
    #drift part
    Bi=matrix(c(delta,0),nrow=2,ncol=1,byrow=TRUE)
    Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
    
    #random part
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi_tilde)
    
    
    #hidden states link matrix
    Ti=matrix(c(1,1/beta*(1-exp(-beta*delta)),0,exp(-beta*delta)),nrow=2,byrow=TRUE)
    Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
    
    alpha=Ti_tilde%*%alpha+Bi_tilde%*%mu+eta
    
    new_pos=1000*c(alpha[1],alpha[3])
    
    if (is_in_land(st_point(new_pos),land)) {
      new_pos=nearest_shore_point(st_point(new_pos),coastline)
      new_speed=(new_pos-prev_pos)/delta
      
      #update hidden state position
      alpha[1]=new_pos[1]/1000
      alpha[3]=new_pos[2]/1000
      
      #update hidden state velocities
      alpha[2]=new_speed[1]/1000-mu1s[i-1]
      alpha[4]=new_speed[2]/1000-mu2s[i-1]
      
    }
    
    if (is.null(log_sigma_obs)) {
      yobs=c(alpha[1],alpha[3])
    }
    else {
      #measurement error
      epsilon=mvrnorm(2,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #state obs link matrix
      Z=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
      Z_tilde=bdiag(replicate(2,Z,simplify=FALSE))
      
      #noisy observation
      yobs=Z_tilde%*%alpha+epsilon
    }
    
    data_sim=rbind(data_sim,c(yobs[1],yobs[2]))
  }
  data_sim$time=times
  
  
  return (data_sim[,c("y1","y2","time")])
}



#' @description Simulate from a 2D constrained CVM :
#' the mean velocity mu is equal to the shoreline vector, defined by the difference between
#' two consecutive nearest shore points : 
#' \[\tau(D_{shore})=delta0+(tau0-delta)/2*(1-tanh(lambda(D_{shore}-D0))) \]
#' @param mu1 first component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param times vector of sampling times
#' @param log_sigma_obs if not NULL, gaussian noise with standard deviation exp(log_sigma_obs) is added to the sampled locations 
#' @param land sf object (list of polygons) defining the land
#' @param verbose if TRUE, display computations along the way.
#' @return dataframe with dim (n,3) where n is the length of 'times' argument and columns "y1", "y2" are the sampled
#' observations (with or without noise) and "time" is the time of the samples.

sim_constrainedCVM=function(delta0,tau0,D0,lambda,nu,v0,x0,times,log_sigma_obs=NULL,land,verbose=TRUE) {
  
  #number of measurements
  n=length(times)
  
  #initialize matrix for nearest shore points
  nearest_points=matrix(rep(NA,2*n),ncol=2)
  
  nus=nu
  
  if (length(nu)==1) {
    nus=rep(nu,n-1)
  }
  
  alpha=matrix(c(x0[1],v0[1],x0[2],v0[2]),nrow=4,ncol=1,byrow=TRUE)
  
  if (is.null(log_sigma_obs)) {
    yobs=c(alpha[1],alpha[3])
  }
  else {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
    
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  
  data_sim=data.frame("y1"=yobs[1],"y2"=yobs[2])
  
  
  if (!(verbose)) {
    #progress bar
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  for (i in 2:n) {
    
    if (!(verbose)) {
      pb$tick()
    }
    
    #time step
    delta=times[i]-times[i-1]
    
    #previous position in metres
    prev_pos=1000*c(alpha[1],alpha[3])
    
    #parameters
    nearest_points[i-1,]=nearest_shore_point(st_point(prev_pos),land)/1000
    n=prev_pos-nearest_points[i-1,]
    Dshore=sqrt(n[1]^2+n[2]^2)
    tau=delta0+(tau0-delta0)/2*(1+tanh(lambda*(Dshore-D0)))
    nu=nus[i-1]
    beta=1/tau
    sigma=2*nu/sqrt(pi*tau)
    
    if (i==2) {
      mu=matrix(c(0,0),ncol=1,nrow=2,byrow=TRUE)
    }
    else {
      mu=nearest_points[i-1,]-nearest_points[i-2,]
    }
    #covariance matrix
    var_xi=sigma^2/beta^2*(delta+2/beta*(exp(-beta*delta)-1)+1/2/beta*(1-exp(-2*beta*delta)))
    var_zeta=sigma^2/2/beta*(1-exp(-2*beta*delta))
    cov=sigma^2/2/beta^2*(1+exp(-2*beta*delta)-2*exp(-beta*delta))
    Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
    Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
    
    #drift part
    Bi=matrix(c(delta,0),nrow=2,ncol=1,byrow=TRUE)
    Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
    
    #random part
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi_tilde)
    
    
    #hidden states link matrix
    Ti=matrix(c(1,1/beta*(1-exp(-beta*delta)),0,exp(-beta*delta)),nrow=2,byrow=TRUE)
    Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
    
    alpha=Ti_tilde%*%alpha+Bi_tilde%*%mu+eta
    
    if (is.null(log_sigma_obs)) {
      yobs=c(alpha[1],alpha[3])
    }
    else {
      #measurement error
      epsilon=mvrnorm(2,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #state obs link matrix
      Z=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
      Z_tilde=bdiag(replicate(2,Z,simplify=FALSE))
      
      #noisy observation
      yobs=Z_tilde%*%alpha+epsilon
    }
    
    data_sim=rbind(data_sim,c(yobs[1],yobs[2]))
  }
  data_sim$time=times
  
  
  return (data_sim[,c("y1","y2","time")])
}




#' @description Simulate from a 2D reflected CVM (as in Hanks2017- reflected SDE)
#' 
#' @param mu1 first component of the drift parameter, scalar 
#' @param mu2 second component of the drift parameter, scalar
#' @param beta autocorrelation parameter, scalar 
#' @param sigma volatility parameter, scalar
#' @param x0 initial velocity
#' @param v0 initial location
#' @param h constant time step between samples
#' @param Tmax maximum time of sample in hours
#' @param land sf object (list of polygons) defining the land
#' @param coastline sf object (list of polygons and lines) defining the coastline or land
#' @param verbose if TRUE, display computations along the way.
#' @return dataframe with dim (n,3) where n is the length of 'times' argument and columns "x1", "x2" are the sampled
#' observations (with or without noise) and "times" is the time of the samples

sim_HANKSreflected=function(mu1,mu2,beta,sigma,x0,v0,h,Tmax,land,coastline,verbose=FALSE) {
  
  N=floor(Tmax/h)
  
  #initialization
  x1=as.numeric(x0+h*v0)
  x=data.frame("x1"=c(x0[1],x1[1]),"x2"=c(x0[2],x1[2]))
  
  #drift
  mu=c(mu1,mu2)
  
  #progress bar
  if (!(verbose)) {
  pb <- progress_bar$new(
    format = "[:bar] :percent",
    total = N, clear = FALSE, width = 60)
  }
  
  for (i in 2:N) {
    
    if (!(verbose)) {
      pb$tick()
    }
    
    xprev1=as.numeric(x[i,])
    xprev2=as.numeric(x[i-1,])
    
    if (verbose) {
      cat("previous position :",xprev1,"\n")
    }
    w=mvrnorm(1,mu =c(0,0), Sigma =h^3*sigma^2*diag(2))
    xtilde=xprev1*(2-beta*h)+xprev2*(beta*h-1)+beta*h*h*mu+w
    
    if ((is_in_land(st_point(1000*c(xtilde)),land))) {
      if (verbose) {
        cat("new position ",xtilde,"is in LAND \n")
      }
      projection=c(nearest_shore_point(st_point(1000*xtilde),coastline)/1000)
      xtilde=projection
    }
    
    x=rbind(x,c(xtilde))
  }
  
  x$time=seq(0,N*h,by=h)
  return (x)
}



#' @description Euler simulation of a 2D constrained CVM where mu is replace by a vector parallell to the shore on each time step
#' \[v(t+\Delta)=exp(-beta \frac{D0}{Dshore} \Delta)\times v(t)+(1-exp(-beta \frac{D0}{Dshore}\times \Delta))\times R \times n\]
#' where $n$ is the unit normal vector between the shore and the last position and $R$ is rotation a matrix that makes n parallel to the shore
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma  volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param D0 reference distance for the effect of shore on the movement
#' @param x0 initial velocity
#' @param v0 initial location
#' @param h constant small time step between samples
#' @param Tmax maximum time of sample in hours
#' @param land sf object (list of polygons) defining the land
#' @param verbose if TRUE, display computations along the way.
#' @return dataframe with dim (n,3) where n is the length of 'times' argument and columns "x1", "x2" are the sampled
#' observations (with or without noise) and "times" is the time of the samples

sim_normalCVM=function(beta,sigma,D0,x0,v0,h,Tmax,land,verbose=FALSE) {
  
  times=seq(0,Tmax,by=h)
  n=length(times)
  
  betas=beta
  sigmas=sigma
  
  if (length(beta)==1) {
    betas=rep(beta,n-1)
  }
  if (length(sigma)==1) {
    sigmas=rep(sigma,n-1)
  }
  
  par=data.frame(beta=betas,sigmas=sigma)
  data_sim=data.frame("x1"=x0[1],"x2"=x0[2],"v1"=v0[1],"v2"=v0[2])
  
  if (!(verbose)) {
    #progress bar
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  for (i in 2:n) {
    
    if (!(verbose)) {
      pb$tick()
    }
  
    
    #time step
    delta=times[i]-times[i-1]
  
    #previous pos and velocity
    xmoins=c(data_sim[i-1,"x1"],data_sim[i-1,"x2"])
    vmoins=c(data_sim[i-1,"v1"],data_sim[i-1,"v2"])
    
    if (verbose) {
      cat("position :",xmoins,"\n")
      cat("velocity",vmoins,"\n")
    }
    
    #parameters
    beta=par$beta[i-1]
    sigma=par$sigma[i-1]
    
    #nearest shore points
    p=c(nearest_shore_point(st_point(1000*xmoins),land)/1000)
    
    n=xmoins-p
    Dshore=sqrt(n[1]^2+n[2]^2)
    
    if (Dshore==0) {
      data_sim$time=times[1:(i-1)]
      return (data_sim[,c("x1","x2","time")])
    }
    
    R=matrix(c(0,-1,1,0),byrow=TRUE,ncol=2)
    epsilon=mvrnorm(1,mu =c(0,0), Sigma =sigma^2*(1-exp(-2*beta*delta))/(2*beta)*diag(2))
    vplus=exp(-beta*D0/Dshore*delta)*vmoins+(1-exp(-beta*D0/Dshore*delta))*R%*%n/Dshore+epsilon
    xplus=xmoins+delta*vmoins
    
    
    
    data_sim=rbind(data_sim,c(xplus[1],xplus[2],vplus[1],vplus[2]))
  }
  
  data_sim$time=times
  
  #return only position
  return (data_sim[,c("x1","x2","time")])
  
}


#' @description Euler simulation of RACVM 
#' 
#' @param mu1 first component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param omega rotational velocity parameter, scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param times vector of sampling times
#' @param verbose if TRUE, display computations along the way.
#' @return datafram with dim (n,3) where n is the length of 'times' argument and columns "x1", "x2" are the sampled
#' locations and "time".

sim_eulerRACVM=function(mu1,mu2,beta,sigma,omega,v0,x0,h,Tmax,verbose=FALSE,keep_true_pos=TRUE) {
  
  times=seq(0,Tmax,by=h)
  n=length(times)
  
  betas=beta
  sigmas=sigma
  omegas=omega
  mu1s=mu1
  mu2s=mu2
  
  if (length(beta)==1) {
    betas=rep(beta,n)
  }
  if (length(sigma)==1) {
    sigmas=rep(sigma,n)
  }
  if (length(omega)==1) {
    omegas=rep(omega,n)
  }
  if (length(mu1)==1) {
    mu1s=rep(mu1,n)
  }
  if (length(mu2)==1) {
    mu2s=rep(mu2,n)
  }
  
  par=data.frame(mu1=mu1s,mu2=mu2s,sigma=sigmas,beta=betas,omega=omegas)
  data_sim=data.frame("x1"=x0[1],"x2"=x0[2],"x1"=x0[1],"x2"=x0[2])
  
  #progress bar
  if (!(verbose)) {
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  for (i in 2:n) {
    
    if (!(verbose)) {
      pb$tick()
    }
    
    #previous pos and velocity
    xmoins=c(data_sim[i-1,"x1"],data_sim[i-1,"x2"])
    vmoins=c(data_sim[i-1,"v1"],data_sim[i-1,"v2"])
    
    if (verbose) {
      cat("velocity :",vmoins,"\n")
      cat("position :",xmoins,"\n")
    }
    
    #parameters
    tau=par$tau[i-1]
    nu=par$nu[i-1]
    omega=par$omega[i-1]
    sigma=2*nu/sqrt(pi*tau)
    mu=matrix(c(par$mu1[i-1],par$mu2[i-1]),ncol=1,nrow=2,byrow=TRUE)
    
    #rotation matrix
    A=matrix(c(1/tau,-omega,omega,1/tau),nrow=2,byrow=TRUE)
    
    #brownian motion increment
    epsilon=mvrnorm(1,mu =c(0,0), Sigma =sigma^2*h*diag(2))
    
    #new velocity and positions
    vplus=vmoins-A%*%(vmoins-mu)*h+epsilon
    xplus=xmoins+vmoins*h
    
    data_sim=rbind(data_sim,c(xplus[1],xplus[2],vplus[1],vplus[2]))
  }
  
  data_sim$time=times
  
  #return only position
  return (data_sim[,c("x1","x2","time")])
}


#' @description Simulate from RACVM with explicit formulas
#' 
#' @param mu1 first component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param mu2 second component of the drift parameter, scalar or vector of scalar with size length(times)-1
#' @param beta autocorrelation parameter, scalar or vector of scalar with size length(times)-1
#' @param sigma volatility parameter, scalar or vector of scalar with size length(times)-1
#' @param omega rotational velocity parameter, scalar or vector of scalar with size length(times)-1
#' @param v0 initial velocity
#' @param x0 initial location
#' @param times vector of sampling times
#' @param log_sigma_obs if not NULL, gaussian noise with standard deviation exp(log_sigma_obs) is added to the sampled locations 
#' @param verbose if TRUE, display computations along the way.
#' @return datafram with dim (n,2) where n is the length of 'times' argument and columns "x1", "x2" are the sampled
#' locations 

sim_RACVM=function(mu1,mu2,beta,sigma,omega,v0,x0,times,log_sigma_obs=NULL,verbose=FALSE,keep_true_pos=FALSE) {

  #number of observations
  n=length(times)
  
  #data frames of parameters values
  betas=beta
  sigmas=sigma
  omegas=omega
  mu1s=mu1
  mu2s=mu2
  
  if (length(beta)==1) {
    betas=rep(beta,n)
  }
  if (length(sigma)==1) {
    sigmas=rep(sigma,n)
  }
  if (length(omega)==1) {
    omegas=rep(omega,n)
  }
  if (length(mu1)==1) {
    mu1s=rep(mu1,n)
  }
  if (length(mu2)==1) {
    mu2s=rep(mu2,n)
  }
  
  par=data.frame(mu1=mu1s,mu2=mu2s,beta=betas,sigma=sigmas,omega=omegas)
  
  #initialization
  
  alpha=matrix(c(x0[1],x0[2],v0[1]-par$mu1[1],v0[2]-par$mu2[1]),nrow=4,ncol=1,byrow=TRUE)
  
  #if log_sigma_obs is not null, add noise
  if (!(is.null(log_sigma_obs))) {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  else {
    yobs=x0
  }
 
  data_sim=data.frame("y1"=yobs[1],"y2"=yobs[2],"x1"=x0[1],"x2"=x0[2])
  
  #progress bar
  if (!(verbose)) {
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
    
  #loop over the observation times
  for (i in 2:n) {
    
    #update progress bar
    if (!(verbose)) {
      pb$tick()
    }
    
    #time step
    delta=times[i]-times[i-1]
    
    #previous pos and velocity
    xmoins=c(data_sim[i-1,"x1"],data_sim[i-1,"x2"])
    vmoins=c(data_sim[i-1,"v1"],data_sim[i-1,"v2"])
    
    if (verbose) {
      cat("position",xmoins,"\n")
      cat("velocity",vmoins,"\n")
    }
    
    #parameters
    mu=matrix(c(par$mu1[i-1],par$mu2[i-1]),ncol=1,nrow=2,byrow=TRUE)
    beta=par$beta[i-1]
    sigma=par$sigma[i-1]
    omega=par$omega[i-1]
  
  
    C=beta^2+omega^2
    A=matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
    invA=1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
    R=matrix(c(cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta)),byrow=TRUE,nrow=2)
    expAdelta=exp(-beta*delta)*R
    
    #covariance matrix
    var_xi=sigma^2/C*(delta+(omega^2-3*beta^2)/(2*beta*C)-exp(-2*delta*beta)/(2*beta)+
                        2*exp(-delta*beta)*(beta*cos(omega*delta)-omega*sin(omega*delta))/C)
    var_zeta=sigma^2/(2*beta)*(1-exp(-2*delta*beta))
    cov1=sigma^2/(2*C)*(1+exp(-2*delta*beta)-2*exp(-delta*beta)*cos(omega*delta))
    cov2=sigma^2/C*(exp(-delta*beta)*sin(omega*delta)-omega/(2*beta)*(1-exp(-2*delta*beta)))
    Qi=matrix(c(var_xi,0,cov1,cov2,0,var_xi,cov2,cov1,cov1,cov2,var_zeta,0,cov2,cov1,0,var_zeta),nrow=4,byrow=TRUE)
    
    
    #random part
    eigenvalues=eigen(Qi)$values
    if (verbose) {
      cat("eigenvalues of covariance matrix :",eigenvalues,"\n")
    }
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi,tol=1e-1)
    
    
    
    #hidden states link matrix
    B1=diag(2)
    B2=invA%*%(diag(2)-expAdelta)
    B3=matrix(c(0,0,0,0),nrow=2)
    B4=expAdelta
    # Combine the matrices into a 4x4 matrix
    Ti <- matrix(0, nrow = 4, ncol = 4)
    Ti[1:2, 1:2] <- B1
    Ti[1:2, 3:4] <- B2
    Ti[3:4, 1:2] <- B3
    Ti[3:4, 3:4] <- B4
    
    if (verbose) {
      cat("Ti :")
      print(Ti)
      cat("\n","eta :")
      print(eta)
    }
    
    #state obs link matrix
    Z=matrix(c(1,0,0,0,0,1,0,0),nrow=2,ncol=4,byrow=TRUE)
    
    #drift part
    Bi<- matrix(0, nrow = 4, ncol = 2)
    Bi[1:2, 1:2] <- delta*diag(2)
    Bi[3:4,1:2] <- 0*diag(2)
    
    
    
    alpha=Ti%*%alpha+Bi%*%mu+eta
    
    if (!(is.null(log_sigma_obs))) {
      
      #measurement error
      epsilon=mvrnorm(1,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #noisy observation
      yobs=Z%*%alpha+epsilon
    }
    else {
      yobs=c(alpha[1],alpha[2])
    }
    
    data_sim=rbind(data_sim,c(yobs[1],yobs[2],alpha[1],alpha[2]))
    
    if (verbose) {
      cat("\n")
    }
  }
  
  data_sim$time=times
  
  if (keep_true_pos) {
    return (data_sim)}
  else {
    #return only position
    return (data_sim[,c("y1","y2","time")])
  }
}


#' @description Simulation of a 2D (simplified) CRCVM with parametrization 

#' @param ftau : smooth function with two arguments Dshore and phi, that gives the value of tau. We can choose it constant.
#' @param fomega : smooth function with two arguments Dshore and phi, that gives the value of tau
#' @param nu volatility parameter: vector with the same size as times, or scalar.
#' @param log_sigma_obs log of standard deviation for measurement error.
#'  Default : NULL (no error)
#' @param x0 initial velocity
#' @param v0 initial location
#' @param times numeric vector of sampling times
#' @param land sf object (list of polygons) defining the land
#' @param verbose if TRUE, display computations along the way.
#' @return list with two components: 
#' - "sim" is dataframe with dim (n,3) where n is the length of 'times' argument and columns "y1", "y2" are the sampled
#' observations (with or without noise) and "times" is the time of the samples. 
#' - shore is a dataframe with columns "p1","p2" that are the coordinates of the nearest points to the shore (useful for the plots)


sim_simplified_CRCVM=function(ftau,fomega,nu,log_sigma_obs=NULL,v0,x0,times,land,verbose=FALSE) {
  
  #number of samples
  n=length(times)
  
  #volatility
  nus=nu
  if (length(nu)==1) {
    nus=rep(nu,n)
  }
  
  
  #initialization
  alpha=matrix(c(x0[1],x0[2],v0[1],v0[2]),nrow=4,ncol=1,byrow=TRUE)
  
  #if log_sigma_obs is not null, add noise
  if (!(is.null(log_sigma_obs))) {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  else {
    yobs=x0
  }
  
  data_sim=data.frame("x1"=x0[1],"x2"=x0[2],"v1"=v0[1],"v2"=v0[2],"y1"=yobs[1],"y2"=yobs[2])
  data_shore=data.frame("p1"=NA,"p2"=NA)
  
  #progress bar
  if (!(verbose)) {
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  #loop over observation times
  for (i in 2:n) {
    
    #update progress bar
    if (!(verbose)) {
      pb$tick()
    }
    
    
    #time step
    delta=times[i]-times[i-1]
    
    #previous pos and velocity
    xmoins=c(data_sim[i-1,"x1"],data_sim[i-1,"x2"])
    vmoins=c(data_sim[i-1,"v1"],data_sim[i-1,"v2"])
    
    if (verbose) {
      cat("new position",xmoins,"\n")
      cat("new velocity",vmoins,"\n")
    }
    
    #nearest shore point
    p=c(nearest_shore_point(st_point(1000*xmoins),land)/1000)
    
    #normal vector
    normal=xmoins-p
    
    #distance to boundary
    Dshore=sqrt(normal[1]^2+normal[2]^2)
    
    
    #if point is already on the land, stop simulation
    if (Dshore==0) {
      data_sim$time=times[1:(i-1)]
      l=list("sim"=data_sim[1:(i-2),c("y1","y2","time")],"shore"=data_shore[2:(i-1),])
      return (l)
    }
    
    # Calculate the signed angle in radians.
    theta=signed_angle(normal,vmoins)
    
    #deviation angle
    if (0<=theta && theta <=pi/2) {
      phi=pi/2-theta
    }
    else if (-pi/2<theta && theta< 0) {
      phi=-(theta+pi/2)
    }
    else {
      phi=0
    }
    
    #add to shore dataframe
    data_shore=rbind(data_shore,c(as.numeric(p)))
    
    #parameters
    tau=ftau(Dshore,phi)
    omega=fomega(Dshore,phi)
    sigma=sigmas[i-1]

    
    if (verbose) {
      cat("Distance to shore :", Dshore,"\n")
      cat("Angular velocity :",omega,"\n")
      cat("Persistence =",tau,"\n")
    }
    
    C=1/tau^2+omega^2
    A=matrix(c(1/tau,-omega,omega,1/tau),nrow=2,byrow=TRUE)
    invA=1/C*matrix(c(1/tau,omega,-omega,1/tau),nrow=2,byrow=TRUE)
    R=matrix(c(cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta)),byrow=TRUE,nrow=2)
    expAdelta=exp(-delta/tau)*R
    
    #covariance matrix
    var_xi=sigma^2/C*(delta+(omega^2-3/tau^2)/(2/tau*C)-exp(-2*delta/tau)/(2/tau)+
                        2*exp(-delta/tau)*(1/tau*cos(omega*delta)-omega*sin(omega*delta))/C)
    var_zeta=sigma^2*tau/2*(1-exp(-2*delta/tau))
    cov1=sigma^2/(2*C)*(1+exp(-2*delta/tau)-2*exp(-delta/tau)*cos(omega*delta))
    cov2=sigma^2/C*(exp(-delta/tau)*sin(omega*delta)-omega/(2/tau)*(1-exp(-2*delta/tau)))
    Qi=matrix(c(var_xi,0,cov1,cov2,0,var_xi,cov2,cov1,cov1,cov2,var_zeta,0,cov2,cov1,0,var_zeta),nrow=4,byrow=TRUE)
    
    
    #random part
    if (verbose) {
      cat("eigenvalues of covariance matrix :",eigen(Qi)$values,"\n")
    }
    
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi,tol=1e-1)
    
    
    #hidden states link matrix
    B1=diag(2)
    B2=invA%*%(diag(2)-expAdelta)
    B3=matrix(c(0,0,0,0),nrow=2)
    B4=expAdelta
    # Combine the matrices into a 4x4 matrix
    Ti <- matrix(0, nrow = 4, ncol = 4)
    Ti[1:2, 1:2] <- B1
    Ti[1:2, 3:4] <- B2
    Ti[3:4, 1:2] <- B3
    Ti[3:4, 3:4] <- B4
    
    alpha=Ti%*%alpha+eta
    
    #if no noise
    if (is.null(log_sigma_obs)) {
      yobs=c(alpha[1],alpha[2])
    }
    else {
      
      #measurement error
      epsilon=mvrnorm(1,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #state obs link matrix
      Z=matrix(c(1,0,0,0,0,1,0,0),nrow=2,ncol=4,byrow=TRUE)
      
      #noisy observation
      yobs=Z%*%alpha+epsilon
    }
    
    data_sim=rbind(data_sim,c(alpha[1],alpha[2],alpha[3],alpha[4],yobs[1],yobs[2]))
  }
  
  data_sim$time=times
  
  l=list("sim"=data_sim[1:(n-1),c("y1","y2","time")],"shore"=data_shore[2:n,])
  return (l)
}





#' @description Simulation of a 2D CRCVM with tau and omega as functions of distance to shore and angle theta
#'
#' @param ftau : smooth function with two arguments Dshore and theta, that gives the value of tau. We can choose it constant.
#' @param fomega : smooth function with two arguments Dshore and theta, that gives the value of tau
#' @param nu volatility parameter: vector with the same size as times, or scalar.
#' @param log_sigma_obs log of standard deviation for measurement error.
#'  Default : NULL (no error)
#' @param x0 initial velocity
#' @param v0 initial location
#' @param times numeric vector of sampling times
#' @param land sf object (list of polygons) defining the land
#' @param verbose if TRUE, display computations along the way.
#' @return list with two components: 
#' - "sim" is dataframe with dim (n,3) where n is the length of 'times' argument and columns "y1", "y2" are the sampled
#' observations (with or without noise) and "times" is the time of the samples. 
#' - shore is a dataframe with columns "p1","p2" that are the coordinates of the nearest points to the shore (useful for the plots)


sim_theta_CRCVM=function(ftau,fomega,fnu,log_sigma_obs=NULL,v0,x0,times,land,verbose=FALSE) {
  
  #number of samples
  n=length(times)
  
  
  
  #initialization
  alpha=matrix(c(x0[1],x0[2],v0[1],v0[2]),nrow=4,ncol=1,byrow=TRUE)
  
  #if log_sigma_obs is not null, add noise
  if (!(is.null(log_sigma_obs))) {
    epsilon=mvrnorm(2,mu =c(0,0), Sigma =diag(exp(log_sigma_obs),2))
    yobs=matrix(c(x0[1]+epsilon[1],x0[2]+epsilon[2]),nrow=2,ncol=1,byrow=TRUE)
  }
  else {
    yobs=x0
  }
  
  data_sim=data.frame("x1"=x0[1],"x2"=x0[2],"v1"=v0[1],"v2"=v0[2],"y1"=yobs[1],"y2"=yobs[2])
  data_shore=data.frame("p1"=NA,"p2"=NA,"Dshore"=NA,"theta"=NA)
  
  #progress bar
  if (!(verbose)) {
    pb <- progress_bar$new(
      format = "[:bar] :percent",
      total = n, clear = FALSE, width = 60)
  }
  
  #loop over observation times
  for (i in 2:n) {
    
    #update progress bar
    if (!(verbose)) {
      pb$tick()
    }
    
    
    #time step
    delta=times[i]-times[i-1]
    
    #previous pos and velocity
    xmoins=c(data_sim[i-1,"x1"],data_sim[i-1,"x2"])
    vmoins=c(data_sim[i-1,"v1"],data_sim[i-1,"v2"])
    
    if (verbose) {
      cat("new position",xmoins,"\n")
      cat("new velocity",vmoins,"\n")
    }
    
    #nearest shore point
    p=c(nearest_shore_point(st_point(xmoins),land))
    
    #normal vector
    normal=xmoins-p
    
    #distance to boundary
    Dshore=sqrt(normal[1]^2+normal[2]^2)
    
    
    #if point is already on the land, stop simulation
    if (Dshore==0) {
      data_sim$time=times[1:(i-1)]
      l=list("sim"=data_sim[1:(i-2),c("y1","y2","time")],"shore"=data_shore[2:(i-1),])
      cat("Stopped : process reached land ! \n")
      return (l)
    }
    
    # Calculate the signed angle in radians.
    theta=signed_angle(normal,vmoins)
    
    #add to shore dataframe
    data_shore=rbind(data_shore,c(as.numeric(p),Dshore,theta))
    
    #parameters
    cov_data=data.frame(DistanceShore=Dshore,theta=theta)
    tau=ftau(cov_data)
    omega=fomega(cov_data)
    nu=fnu(cov_data)
    sigma=2*nu/sqrt(pi*tau)
    
    if (verbose) {
      cat("Distance to shore :", Dshore,"\n")
      cat("Angular velocity :",omega,"\n")
      cat("Persistence =",tau,"\n")
    }
    
    C=1/tau^2+omega^2
    A=matrix(c(1/tau,-omega,omega,1/tau),nrow=2,byrow=TRUE)
    invA=1/C*matrix(c(1/tau,omega,-omega,1/tau),nrow=2,byrow=TRUE)
    R=matrix(c(cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta)),byrow=TRUE,nrow=2)
    expAdelta=exp(-delta/tau)*R
    
    #covariance matrix
    var_xi=sigma^2/C*(delta+(omega^2-3/tau^2)/(2/tau*C)-exp(-2*delta/tau)/(2/tau)+
                        2*exp(-delta/tau)*(1/tau*cos(omega*delta)-omega*sin(omega*delta))/C)
    var_zeta=sigma^2*tau/2*(1-exp(-2*delta/tau))
    cov1=sigma^2/(2*C)*(1+exp(-2*delta/tau)-2*exp(-delta/tau)*cos(omega*delta))
    cov2=sigma^2/C*(exp(-delta/tau)*sin(omega*delta)-omega/(2/tau)*(1-exp(-2*delta/tau)))
    Qi=matrix(c(var_xi,0,cov1,cov2,0,var_xi,cov2,cov1,cov1,cov2,var_zeta,0,cov2,cov1,0,var_zeta),nrow=4,byrow=TRUE)
    
    
    #random part
    if (verbose) {
      cat("eigenvalues of covariance matrix :",eigen(Qi)$values,"\n")
    }
    
    eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi,tol=1e-1)
    
    
    #hidden states link matrix
    B1=diag(2)
    B2=invA%*%(diag(2)-expAdelta)
    B3=matrix(c(0,0,0,0),nrow=2)
    B4=expAdelta
    # Combine the matrices into a 4x4 matrix
    Ti <- matrix(0, nrow = 4, ncol = 4)
    Ti[1:2, 1:2] <- B1
    Ti[1:2, 3:4] <- B2
    Ti[3:4, 1:2] <- B3
    Ti[3:4, 3:4] <- B4
    
    alpha=Ti%*%alpha+eta
    
    #if no noise
    if (is.null(log_sigma_obs)) {
      yobs=c(alpha[1],alpha[2])
    }
    else {
      
      #measurement error
      epsilon=mvrnorm(1,mu=c(0,0), Sigma =diag(exp(log_sigma_obs)^2,2))
      
      #state obs link matrix
      Z=matrix(c(1,0,0,0,0,1,0,0),nrow=2,ncol=4,byrow=TRUE)
      
      #noisy observation
      yobs=Z%*%alpha+epsilon
    }
    
    data_sim=rbind(data_sim,c(alpha[1],alpha[2],alpha[3],alpha[4],yobs[1],yobs[2]))
  }
  
  data_sim$time=times
  
  l=list("sim"=data_sim[1:(n-1),c("y1","y2","time")],"shore"=data_shore[2:n,])
  return (l)
}



