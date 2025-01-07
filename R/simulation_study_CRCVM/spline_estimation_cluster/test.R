library(ggplot2)


domain="fjords"
TMAX=12
N_ID=6
DMIN=1
cov="DistanceShore"
pattern<- sprintf("result_%s_%dh_%dID_%dkm_%s[1-9][0-9]*",domain,TMAX,N_ID,DMIN,cov)
files<- list.files(path=domain,pattern = pattern)

#get the set up
set_up_file=paste("set_up_","fjords",".R",sep="")
source(set_up_file)    

angle1=pi/2
angle2=-pi
# Generate the x-values
x_vals <- seq(from=1.5, to=3, length.out=100)

# True functions
omega_away <- fomega_splines(cov_data=data.frame("theta"=rep(angle1, 100), "DistanceShore"=x_vals))
omega_toward <- fomega_splines(cov_data=data.frame("theta"=rep(angle2, 100), "DistanceShore"=x_vals))

# Create initial plots with true functions
plot_omega_away <- ggplot() +
  geom_line(aes(x=x_vals, y=omega_away), color="red") +
  ylab(expression(omega)) +
  xlab(expression(Theta))

plot_omega_toward <- ggplot() +
  geom_line(aes(x=x_vals, y=omega_toward), color="red") +
  ylab(expression(omega)) +
  xlab(expression(Theta))

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
print(plot_omega_away)
print(plot_omega_toward)