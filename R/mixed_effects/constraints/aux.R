
# NECESSARY PACKAGES

#for sde fit
library(smoothSDE)

#for plots
library(ggplot2)
library(plotly)

#to define and extract string formulas
library(splines)
library(formula.tools)
library(stringr)

#for gam models
library(mgcv)

#to save 3D files
library(htmlwidgets)

#for nearest point on coastline
library(sf)

#for block matrix in SSM for simulation
library(Matrix)
library(MASS)


#####################################  FUNCTION TO GET ALL INTERESTING PLOTS OF THE PARAMETERS OF A MODEL  ##########################################################



get_par_plots=function(baseline=NULL,sde,model_name,links=list(),xmins=list(),xmaxs=list(),xlabels=list(),npost=1000,level=0.95) {
  
  
  ####################################################################################################
  #get all the plots of the estimates of the parameters according to the model formulas
  
  #we suppose that each parameter has at least one random effect based on ID covariate
  
  #for each parameter :
  #   - if there are no fixed effects (only random effect),
  #     we only plot the value of the parameter for each narwhal
  #     with its confidence interval
  #   - If there is only one spline fixed effect, we plot the 
  #     fixed effect estimates and the mixed effects on two different plots
  #   - If there are two covariates, two cases are possible :
  #       . either the covariates are orthogonal, in which case we plot the fixed effects and 
  #         the mixed effects for each covariate while forcing the other to 0
  #       . or they are not orthogonal, in which case we plot the fixed effects in 3D
  
  
  #ARGUMENTS :
  #   - baseline is a baseline model to be compared to (default is NULL). 
  #     It Should be given as an entry if we want to plot baseline parameters on the smooth functions graph
  #   - sde is a fitted sde model 
  #   - model_name is a string to put in the name of the saved plots
  #   - links is a list of functions that links the fixed effect covariate and the quantity appearing in the x axis of the plots
  #     The names of the list is a subset of the names of all the covariates appearing in the formulas.
  #   - xmin, xmax : lists of boundary values for each fixed effects covariate in the model to plot the parameters. 
  #      If NULL, then we take the min and max of the observed values for each covariate.
  #   - xlabels is a list of label for the x-axis for each covariate
  #   - npoints is the number of points where the fixed effects are plotted
  
  ################################################################################################
  
  #get data
  data=sde$data()
  
  # get formulas and parameters
  formulas=sde$formulas()
  pars=names(formulas)
  
  #draw coeff from posterior distribution once and for all 
  # to get confidence intervals later
  postcoeff=sde$post_coeff(npost)
  all_coeff_re_post=postcoeff$coeff_re
  all_coeff_fe_post=postcoeff$coeff_fe
  
  #get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  
  #Set link function to identity if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(links))) {
      link=(\(x) x)
      links[[var]]=link
    }
  }
  
  #Set xmins and xmaxs to min and max of observed covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xmins))) {
      xmin=min(data[,var])
      xmins[[var]]=xmin
    }
    if (!(var %in% names(xmaxs))) {
      xmax=max(data[,var])
      xmaxs[[var]]=xmax
    }
  }
  
  #Set xlabels to name of covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xlabels))) {
      xlabels[[var]]=var
    }
  }
  
 
  
  #loop over the formulas of each parameter 
  for (j in seq_along(formulas)) {
    
    #formula for the parameter j
    form=formulas[[j]] 
    
    #name of parameter j
    par=pars[j] 
    
    #fixed effects covariates/variables in the formula
    vars=get_variables(form) 
    vars=vars[vars!="ID"]
    
    #number of covariates
    ncovs=length(vars) 
    
    #if there is no fixed effect
    if (ncovs==0) {
      plot_re_par(sde=sde,par=par,model_name=model_name,npost=npost,level=level,save=TRUE)
    }
    
    #if there are fixed effects
    else {
      
      #if there is only one fixed effect covariable
      if (ncovs==1){
      
      link=links[[vars[1]]]
      xlabel=xlabels[vars[1]]
      xmax=xmaxs[[vars[1]]]
      xmin=xmins[[vars[1]]]
    
      
      plot_fe_par_2D(baseline,sde,model_name,par,npoints=200,xmin,xmax,link,xlabel,
                     all_coeff_re_post,all_coeff_fe_post,level,save=TRUE)
      plot_me_par_2D(baseline,sde,model_name,par,npoints=200,xmin,xmax,link,xlabel,npost,level,save=TRUE)
      }
    
      # if there are two fixed effects covariables
      else if (ncovs==2){
        
        #if the covariates are orthogonal
        if (are_orthogonals(data,vars[1],vars[2])) {
          plot_fe_par_2D_ortho(baseline,sde,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                               all_coeff_re_post,all_coeff_fe_post,level,save=TRUE)
          
          plot_me_par_2D_ortho(baseline,sde,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                               npost,level)
        }
        
      #else, plot in 3D
      else {
        plot_fe_par_3D(baseline,sde,model_name,par,npoints=200,xmins=xmins,xmaxs=xmaxs,links=links,
                       xlabels=xlabels,all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level=0.95,save=TRUE)
        }
      }
    }
  }
}








#####################################################   AUXILIARY FUNCTIONS FOR get_par_plots   ############################################################################


are_orthogonals=function(data,C1,C2){
  inner_product=data[,C1]*data[,C2]
  return (all(inner_product==0))
}




plot_re_par=function(sde,model_name,par,npost=1000,level=0.95,save=TRUE) {
  
  
  ############################################################################
  
  #plot the estimated parameter par for each individual with confidence intervals
  
  #sde is a fitted sde model
  
  #par is a string matching with one parameter name in sde 
  #THAT HAS ONLY RANDOM EFFECT IN THE FORMULA
  
  #############################################################################
  
  #get data and IDS
  data=sde$data()
  IDs=unique(data$ID)
  n_id=length(IDs)
  
  #define ID covariates to get the parameters estimate from
  new_data=data[1:n_id,]
  new_data$ID=IDs
  
  #all parameters estimates
  sde_par <- sde$par(t = "all",new_data=new_data)
  
  #all parameters CIs
  sde_CI <- sde$CI_pointwise(t = "all",new_data=new_data,n_post=npost,level=level)
  
  # Data frame for point estimates
  sde_par_df <- data.frame(ID=new_data$ID,par =sde_par[, par])
  
  # Data frame for CIs
  sde_ci_df <- data.frame(ID=new_data$ID,lowpar = sde_CI[par, "low",],
                            uppar = sde_CI[par, "upp",])
  #labels for the plot
  labels=paste("A",1:n_id,sep="")
  
  #ggplot
  plotpar=ggplot()+geom_point(data=sde_par_df, aes(ID, par)) +
    geom_errorbar(data=sde_ci_df,aes(x=ID,ymax = uppar, ymin= lowpar), width= 0.1)+
    scale_x_discrete(labels=labels)+
    ylab(par)+
    theme(axis.title.y=element_text(angle=0),axis.text=element_text(size=12),
          axis.title=element_text(size=16,face="bold"))
  
  if (save) {
    if (!(dir.exists(model_name))) {
      dir.create(model_name)
    }
    ggsave(paste(paste("re",model_name,par,sep="_"),".png",sep=""),plot=plotpar,width=10,height=5,path=model_name)
  }
  return (plotpar)
}












plot_fe_par_2D=function(baseline=NULL,sde,model_name,par,npoints=200,xmin=NULL,xmax=NULL,
                        link=(\(x) x),xlabel="x",all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level=0.95,save=TRUE) {
  
  ###########################################################################################
  # plot the fixed effect for a parameter  that depends only on one covariate 
  # through splines in a sde model
  
  #ARGUMENTS
  
  # baseline : a fitted baseline sde model without fixed effect. If not NULL, the baseline values
  #            appear on the plots
  #
  # sde : a fitted (response) sde 
  #
  # model_name : a string to put in the name of the saved plot
  
  # par : the parameter we want to plot
  #
  # npoints : number of points for the plot
  #
  # xmnin : minimum value of the covariate to plot (if null, take minimum of observed values)
  #
  # xmaxs : maximum value of the covariate to plot (if null, take maximum of observed values)
  #
  # link: function to link the covariate to the quantity we want to have on the x-axis in the plots
  #
  # xlabel : label for the quantity in the xaxis for the plot
  #
  # all_coeff_re_post, all_coeff_fe_post: dataframe of samples of the estimated 
  #         re coeffs (as given in sde$post_coeff()). If null, CIs are not plotted.
  #
  # level: level for confidence intervals
  #
  # save : to save or not to save the plots
  
  
  
  ###########################################################################################
  
  #get data
  data=sde$data()
  
  #number of posteriori draws for CIs
  npost=length(all_coeff_re_post[,1])
  
  #get the formula for the parameter of interest
  formulas=sde$formulas()
  form=formulas[[par]]
  
  #number of parameters 
  n_par=length(names(formulas))
  
  #get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  #get index of the parameter of interest in the vector of all parameters
  j=which(names(formulas)==par)[[1]]
  
  #list of strings for each term in the formula for the parameter of interest
  terms=colnames(attr(terms(form),which="factors")) 
  
  #fixed effects covariates/variables in the formula for the parameter of interest
  vars=get_variables(form) 
  vars=vars[vars!="ID"]
  
  #there is only one covariate
  var=vars[[1]]
  
  #the only fixed effect term is supposed to be at the beginning in the formula
  term=terms[1]
  
  #initialize data frame with covariate values 
  #if boundary values are not defined, take min and max of observations
  if (is.null(xmin)){
    xmin=min(data[,var])
  }
  if (is.null(xmax)){
    xmax=max(data[,var])
  }
  
  #define the covariates values 
  new_data<- data.frame(matrix(0, nrow = npoints, ncol = length(all_vars)))
  # Assign column names to the data frame
  colnames(new_data) <- all_vars
  
  
  #Asssign meaningful values to the covariate that appear in the formula 
  #for the parameter of interest
  new_data[,var]=seq(from=xmin,to=xmax,length.out=npoints)
  
  #Plug in any value for ID since we don't consider random effects
  new_data$ID=rep(data$ID[1],npoints)
  

  #extract spline degree of freedom
  str=paste(par,paste("s(",var,")",sep=""),sep=".")
  q<- as.integer(sde$terms()$ncol_re[str])+1
  #q <- as.integer(str_extract(term, "(?<=k\\s=\\s)\\d+"))
  
  #get model matrices 
  mats=sde$make_mat(new_data=new_data)
  X_fe=mats$X_fe
  X_re=mats$X_re
  
  #extract fixed effects coeff names
  fe_names=paste(rep(paste(par,".s(",var,").",sep=""),q-1),1:(q-1),sep="")
  
  #keep coeffs and columns matching those names in the random effects
  coeff_re_names=rownames(sde$coeff_re())
  index=which(coeff_re_names %in% fe_names)
  X_re=X_re[,index]
  coeff_fe=sde$coeff_fe()
  coeff_re=sde$coeff_re()[index,]
  
  #linear predictor
  lp=X_fe%*%coeff_fe+X_re%*%coeff_re
  lp=matrix(lp,ncol=n_par)
  
  #parameters on response scale
  par_mat <- matrix(sapply(1:ncol(lp), function(i) {
    sde$invlink()[[i]](lp[,i])
  }), ncol = ncol(lp))
  colnames(par_mat) <- names(sde$invlink())
  
  #dataframe for fixed effects estimations
  est=as.data.frame(par_mat[,par])
  est$cov=new_data[,var]
  est$X=link(new_data[,var])
  colnames(est)=c("par",var,"X")
  
  #if there is a baseline, add its values to the df
  if (!(is.null(baseline))) {
    
    #get intercept value on parameter scale
    coeff_fe_baseline=baseline$coeff_fe()[j,1]
    link_fn=sde$invlink()[[j]]
    intercept=link_fn(coeff_fe_baseline)
    est$par_baseline=rep(intercept,npoints)
    
    #draws from posterior distribution of the intercepts
    sd_coeff_fe_baseline=as.list(baseline$tmb_rep(),what="Std")$coeff_fe[j,1]
    if (is.nan(sd_coeff_fe_baseline)) {
      sd_coeff_fe_baseline=0
    }
    par_baseline_post=link_fn(rnorm(npost,mean=coeff_fe_baseline,sd=sd_coeff_fe_baseline))
    alpha <- (1 - level)/2
    quantiles=quantile(par_baseline_post,probs = c(alpha, 1 - alpha))
    
    est$lowpar_baseline=rep(quantiles[1],npoints)
    est$uppar_baseline=rep(quantiles[2],npoints)
    
  }
  
  #differential of population smooth
  h=est$X[2]-est$X[1]
  X_diff=est$X[2:(npoints-1)]
  par_diff=(est$par[1:(npoints-2)]-est$par[3:npoints])/(2*h)
  est_diff=est[2:(npoints-1),]
  est_diff$par=par_diff
  est_diff$X=X_diff
  pdiff=ggplot()+geom_line(data=est_diff,aes(X,par),col="black")+xlab(xlabel)+ylab(paste("derivative of",par))
  
  
  #if we are given posterior draws of the parameters, add CIs to the df 
  if (!(is.null(all_coeff_fe_post)) & !(is.null(all_coeff_re_post))) {
  
    #confidence intervals for the population smooth
    coeff_re_post=all_coeff_re_post[,index]
    coeff_fe_post=all_coeff_fe_post
    
    
    n = nrow(X_fe)/n_par
   
    
    #Get corresponding SDE parameters
    post_par <- array(sapply(1:npost, function(i) {
      sde$par(t = "all", 
                X_fe = X_fe, X_re = X_re, 
                coeff_fe =coeff_fe_post[i,],
                coeff_re = coeff_re_post[i,],
                resp = TRUE)  
    }), dim = c(n, n_par, length(all_coeff_re_post[,1])))
    
    dimnames(post_par)[[2]] <- names(sde$invlink())
    
    alpha <- (1 - level)/2
    sde_CI <- apply(post_par, c(1, 2), quantile, probs = c(alpha, 1 - alpha))
    
    sde_CI <- aperm(sde_CI, c(3, 1, 2))
    dimnames(sde_CI)[[2]] <- c("low", "upp")
    
    # Add 95% quantiles of posterior draws to the dataframe
    est$lowpar <- sde_CI[par, "low",]
    est$uppar<-sde_CI[par, "upp",]
    
  }
  
  #only parameter estimates
  if (is.null(baseline) & (is.null(all_coeff_fe_post)) & (is.null(all_coeff_re_post))) {
  
    p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabel)+ylab(par)
  }
  
  #parameter estimates and baseline
  if (!(is.null(baseline)) & (is.null(all_coeff_fe_post)) & (is.null(all_coeff_re_post))) {
    
    p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabel)+ylab(par)+
      geom_line(aes(X,par_baseline),linetype="dashed",col="black")+
      geom_ribbon(aes(x=X,ymin=lowpar_baseline,ymax=uppar_baseline),fill="grey",alpha=0.4)
  }
  
  #parameter estimates and CIs
  if (is.null(baseline) & !(is.null(all_coeff_fe_post)) & !(is.null(all_coeff_re_post))) {
    
    p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabel)+ylab(par)+
      geom_ribbon(aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)
    
  }
  
  #parameter estimates, baseline and CIs
  else {
    p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabel)+ylab(par)+
      geom_line(aes(X,par_baseline),linetype="dashed",col="black")+
      geom_ribbon(aes(x=X,ymin=lowpar_baseline,ymax=uppar_baseline),fill="grey",alpha=0.4)+
      geom_ribbon(aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)
    
  }
  
  #save plots
  if (save) {
    if (!(dir.exists(model_name))) {
      dir.create(model_name)
    }
    ggsave(paste(paste("fe",model_name,par,var,sep="_"),".png",sep=""),plot=p,width=10,height=5,path=model_name)
    ggsave(paste(paste("diff_fe",model_name,par,var,sep="_"),".png",sep=""),plot=pdiff,width=10,height=5,path=model_name)
  }
  
  return (p)
}











plot_me_par_2D=function(baseline=NULL,sde,model_name,par,npoints=200,xmin=NULL,xmax=NULL,link=(\(x) x),xlabel="x",npost=1000,level=0.95,save=TRUE) {
  
  ########################################################################################################################
  
  #function to plot the estimates for "par" for each individual (taking into account random effects) 
  # when there is only one fixed effect covariate
  
  #ARGUMENTS
  
  #baseline is a baseline model without fixed effects to be compared to (default is NULL). 
  
  #sde is a fitted sde model 
  
  #par is the name of the parameter we want to plot
  
  #model_name is a string to put in the name of the saved plots
  
  #link is a function that links the covariate and the quantity appearing in the x axis of the plots
  
  #xmin, xmax : boundary values of the covariates to plot the parameters. If NULL, then we take the min and max
  #of the observed values of each covariate
  
  #xlabel is the label of the x-axis in the plots
  
  #npoints is the number of points where the smooth population function is plotted
  
  #save: boolean
  
  
  #######################################################################################################################
  
  
  #get data
  data=sde$data()
  
  #get the formula for the parameter of interest
  formulas=sde$formulas()
  form=formulas[[par]]
  
  #get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  #fixed effects covariates/variables in the formula for the parameter of interest
  vars=get_variables(form) 
  vars=vars[vars!="ID"]
  
  #there is only one covariate
  var=vars[[1]]
  
  #IDs
  IDs=unique(data$ID)
  n_id=length(IDs)
  
  
  #set values for current covariate
  #if boundary values are not defined, take min and max of observations
  if (is.null(xmin)){
    xmin=min(data$var)
  }
  if (is.null(xmax)){
    xmax=max(data$var)
  }
  
  
  #define the covariates values 
  new_data<- data.frame(matrix(0, nrow = npoints*n_id, ncol = length(all_vars)))
  
  # Assign column names to the data frame
  colnames(new_data) <- all_vars
  
  #Asssign meaningful values to the covariate that appear in the formula 
  #for the parameter of interest
  cov_values=seq(from=xmin,to=xmax,length.out=npoints)
  cov_values=cov_values[rep(seq_len(npoints),n_id)]
  new_data[,var]=cov_values
  new_data$ID=IDs[rep(seq_len(n_id), each = npoints)] 
  
  sde_par <- sde$par(t = "all",new_data=new_data)
  sde_CI <- sde$CI_pointwise(t = "all",new_data=new_data,n_post=npost,level=level)
  
  
  # Data frame for point estimates
  sde_par_df <- data.frame(ID=new_data$ID,Cov=new_data[,var],par_estimates = sde_par[, par])
  
  
  # Data frame for CIs
  sde_ci_df <- data.frame(ID=new_data$ID,Cov=new_data[,var],
                            lowpar = sde_CI[par, "low",],
                            uppar= sde_CI[par, "upp",])
  
  
  sde_par_df$X=link(new_data[,var])
  sde_ci_df$X=link(new_data[,var])
  
  #labels for the plot
  labels=paste("A",1:n_id,sep="")
  
  p=ggplot()+geom_line(aes(X, par_estimates,col=ID),data=sde_par_df) +
    geom_ribbon(data=sde_ci_df,aes(x=X,ymin=lowpar,ymax=uppar,fill=ID),alpha=0.2)+
    scale_color_discrete(name="ID",labels=labels)+
    scale_fill_discrete(name="ID",labels=labels)+xlab(xlabel)+ylab(par)
  
  
  if (!(is.null(baseline))) {
    baseline_par=baseline$par(t = "all",new_data=new_data)
    baseline_CI=baseline$CI_pointwise(t="all",new_data=new_data)
    baseline_par_df=data.frame(ID=new_data$ID,Cov=new_data[,var],par_baseline = baseline_par[, par])
    baseline_par_df$X=link(new_data[,var])
    p=p+geom_line(aes(X, par_baseline,col=ID),data=baseline_par_df,linetype="dashed")
  }
  
  if (save) {
    
    if (!(dir.exists(model_name))) {
      dir.create(model_name)
    }
    ggsave(paste(paste("me",model_name,par,var,sep="_"),".png"),plot=p,width=10,height=5,path=model_name)
  }
  
  return (p)
}



plot_fe_par_2D_ortho=function(baseline=NULL,sde,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                              all_coeff_re_post,all_coeff_fe_post,level=0.95,save=TRUE) {
  
  
  #################################################################################################################
  
  # plot the fixed effect for a parameter that depends on two orthogonal covariates through splines
  # a plot of the values of the parameter according to each link(covariate) is produced
  #random effects need to be the last term in the formulas
  
  
  #ARGUMENTS
  # baseline is a fitted sde without covariates to compare to. Intercepts of the baseline are displayed on the plots.
  #
  # sde is a fitted sde with mixed effects
  #
  #model_name is a string that appears in the name of the plot saved
  #
  #link is a function that links the covariate and the quantity appearing in the x axis of the plots
  #
  #xlabel is the label of the x-axis in the plots
  #
  #xmin, xmax boundary values of the covariates to plot the parameters. If NULL, then we take the min and max
  #of the observed values of each covariate
  #
  #npoints is the number of points where the smooth population function is plotted
  #
  #level is the level of the confidence interval
  #
  
  ####################################################################################################################
  
  
  #get data
  data=sde$data()
  
  #number of posterior draws from the parameters
  npost=length(all_coeff_re_post[,1])
  
  
  #get the formula for the parameter of interest
  formulas=sde$formulas()
  form=formulas[[par]]
  
  #number of parameters
  n_par=length(names(formulas))
  
  #get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  #get index of the parameter of interest in the vector of all parameters
  j=which(names(formulas)==par)[[1]]
  
  #list of strings for each term in the formula for the parameter of interest
  terms=colnames(attr(terms(form),which="factors")) 
  
  #fixed effects covariates/variables in the formula for the parameter of interest
  vars=get_variables(form) 
  vars=vars[vars!="ID"]
  
  #initialize data frame with covariate values 
  new_data<- data.frame(matrix(0, nrow = npoints, ncol = length(all_vars)))
  
  # Assign column names to the data frame
  colnames(new_data) <- all_vars
  
  #if boundaries for covariate values are no given in xmins, xmaxs arguments
  #take min and max of the observations
  if (is.null(xmins)) {
    xmins=c(min(data[,vars[1]]),min(data[,vars[2]]))
  }
  if (is.null(xmaxs)) {
    xmaxs=c(max(data[,vars[1]]),max(data[,vars[2]]))
  }
  
  list_plots=list()
  
  #loop over the two orthogonal covariates
  for (i in 1:2) {
    
    var=vars[i]
    term=terms[i]
    
    #extract spline degree of freedom
    str=paste(par,paste("s(",var,")",sep=""),sep=".")
    q<- as.integer(sde$terms()$ncol_re[str])+1
    
    #q <- as.integer(str_extract(term, "(?<=k\\s=\\s)\\d+"))
    
    
    new_data[,var]=seq(from=xmins[[var]],to=xmaxs[[var]],length.out=npoints)
    
    #get model matrices 
    mats=sde$make_mat(new_data=new_data)
    X_fe=mats$X_fe
    X_re=mats$X_re
    
    #extract fixed effects coeff names
    fe_names=paste(rep(paste(par,".s(",var,").",sep=""),q-1),1:(q-1),sep="")
    
    #keep coeffs and columns matching those names in the random effects
    coeff_re_names=rownames(sde$coeff_re())
    index=which(coeff_re_names %in% fe_names)
    X_re=X_re[,index]
    coeff_fe=sde$coeff_fe()
    coeff_re=sde$coeff_re()[index,]
    
    #linear predictor
    lp=X_fe%*%coeff_fe+X_re%*%coeff_re
    lp=matrix(lp,ncol=n_par)
    
    #parameters on response scale
    par_mat <- matrix(sapply(1:ncol(lp), function(i) {
      sde$invlink()[[i]](lp[,i])
    }), ncol = ncol(lp))
    colnames(par_mat) <- names(sde$invlink())
    
    #dataframe for fixed effects estimations
    est=as.data.frame(par_mat[,par])
    est$cov=new_data[,i]
    est$X=links[[var]](new_data[,var])
    colnames(est)=c("par","cov","X")
    
    #add baseline intercept value to the plot
    if (!(is.null(baseline))) {
      
      #get intercept value on parameter scale
      coeff_fe_baseline=baseline$coeff_fe()[j,1]
      link_fn=sde$invlink()[[j]]
      intercept=link_fn(coeff_fe_baseline)
      est$par_baseline=rep(intercept,npoints)
      
      #draws from posterior distribution of the intercepts
      sd_coeff_fe_baseline=as.list(baseline$tmb_rep(),what="Std")$coeff_fe[j,1]
      if (is.nan(sd_coeff_fe_baseline)) {
        sd_coeff_fe_baseline=0
      }
      par_baseline_post=link_fn(rnorm(npost,mean=coeff_fe_baseline,sd=sd_coeff_fe_baseline))
      alpha <- (1 - level)/2
      quantiles=quantile(par_baseline_post,probs = c(alpha, 1 - alpha))
      
      est$lowpar_baseline=rep(quantiles[1],npoints)
      est$uppar_baseline=rep(quantiles[2],npoints)
    }
    
    #differential of population smooth
    h=est$X[2]-est$X[1]
    X_diff=est$X[2:(npoints-1)]
    par_diff=(est$par[1:(npoints-2)]-est$par[3:npoints])/(2*h)
    est_diff=est[2:(npoints-1),]
    est_diff$par=par_diff
    est_diff$X=X_diff
    pdiff=ggplot()+geom_line(data=est_diff,aes(X,par),col="black")+xlab(xlabels[[var]])+ylab(paste("derivative of",par))
    
    
    #if we are given posterior draws of the parameters, we can plot the CIs
    if (!(is.null(all_coeff_fe_post)) & !(is.null(all_coeff_re_post))) {
      
      #confidence intervals for the population smooth
      coeff_re_post=all_coeff_re_post[,index]
      coeff_fe_post=all_coeff_fe_post
      
      n = nrow(X_fe)/n_par
      
      #Get corresponding SDE parameters
      post_par <- array(sapply(1:npost, function(i) {
        sde$par(t = "all", 
                  X_fe = X_fe, X_re = X_re, 
                  coeff_fe =coeff_fe_post[i,],
                  coeff_re = coeff_re_post[i,],
                  resp = TRUE)  
      }), dim = c(n, n_par, npost))
      
      dimnames(post_par)[[2]] <- names(sde$invlink())
      
      alpha <- (1 - level)/2
      sde_CI <- apply(post_par, c(1, 2), quantile, probs = c(alpha, 1 - alpha))
      
      sde_CI <- aperm(sde_CI, c(3, 1, 2))
      dimnames(sde_CI)[[2]] <- c("low", "upp")
      
      # Add quantiles of the posterior samples to the df
      est$lowpar=sde_CI[par, "low",]
      est$uppar=sde_CI[par, "upp",]
      
    }
    
    
    #only parameter estimates
    if (is.null(baseline) & (is.null(all_coeff_fe_post)) & (is.null(all_coeff_re_post))) {
      
      p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabels[[var]])+ylab(par)
    }
    
    #parameter estimates and baseline
    if (!(is.null(baseline)) & (is.null(all_coeff_fe_post)) & (is.null(all_coeff_re_post))) {
      
      p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabels[[var]])+ylab(par)+
        geom_line(aes(X,par_baseline),linetype="dashed",col="black")+
        geom_ribbon(aes(x=X,ymin=lowpar_baseline,ymax=uppar_baseline),fill="grey",alpha=0.4)
    }
    
    #parameter estimates and CIs
    if (is.null(baseline) & !(is.null(all_coeff_fe_post)) & !(is.null(all_coeff_re_post))) {
      
      p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabels[[var]])+ylab(par)+
        geom_ribbon(aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)
      
    }
    
    #parameter estimates, baseline and CIs
    else {
      
      p=ggplot(data=est)+geom_line(aes(X,par),col="red")+xlab(xlabels[[var]])+ylab(par)+
        geom_line(aes(X,par_baseline),linetype="dashed",col="black")+
        geom_ribbon(aes(x=X,ymin=lowpar_baseline,ymax=uppar_baseline),fill="grey",alpha=0.4)+
        geom_ribbon(aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)
    }
    
    #save plots
    if (save) {
      
      if (!(dir.exists(model_name))) {
        dir.create(model_name)
      }
      ggsave(paste(paste("fe",model_name,par,var,sep="_"),".png",sep=""),plot=p,width=10,height=5,path=model_name)
      ggsave(paste(paste("diff_fe",model_name,par,var,sep="_"),".png",sep=""),plot=pdiff,width=10,height=5,path=model_name)
    }
    
    list_plots[[var]]=p
  }
  return (list_plots)
}




plot_me_par_2D_ortho=function(baseline=NULL,sde,model_name,par,npoints=200,xmins=list(),xmax=list(),
                              links=list(),xlabels=list(),npost=1000,level=0.95,save=TRUE) {
  
  ##############################################################################################################
  #plot mixed effects of a parameter when there are two orthogonal covariates
  
  ######################################
  
  
  #get data
  data=sde$data()
  
  
  #get the formula for the parameter of interest
  formulas=sde$formulas()
  form=formulas[[par]]
  n_par=length(names(formulas))
  
  ##get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  #fixed effects covariates/variables in the formula for the parameter of interest
  vars=get_variables(form) 
  vars=vars[vars!="ID"]
  
  #IDs
  IDs=unique(data$ID)
  n_id=length(IDs)
  
  
  #initilize dataframe with covariates values 
  new_data<- data.frame(matrix(0, nrow = npoints*n_id, ncol = length(all_vars)))
  
  # Assign column names to the data frame
  colnames(new_data) <- all_vars
  
  #loop over the covariates
  for (i in 1:2) {
    
    #name of the covariate
    var=vars[i]
    
    #Asssign meaningful values to the covariate i that appear in the formula 
    #for the parameter of interest
    cov_values=seq(from=xmins[[var]],to=xmaxs[[var]],length.out=npoints)
    cov_values=cov_values[rep(seq_len(npoints),n_id)]
    new_data[,var]=cov_values
    new_data$ID=IDs[rep(seq_len(n_id), each = npoints)] 
    
    
    sde_par <- sde$par(t = "all",new_data=new_data)
    sde_CI <- sde$CI_pointwise(t = "all",new_data=new_data,n_post=npost,level=level)
    
    
    # Data frame for point estimates
    sde_par_df <- data.frame(ID=new_data$ID,Cov=new_data[,var],par_estimates = sde_par[, par])
    
    
    # Data frame for CIs
    sde_ci_df <- data.frame(ID=new_data$ID,Cov=new_data[,var],
                              lowpar = sde_CI[par, "low",],
                              uppar= sde_CI[par, "upp",])
    
    
    sde_par_df$X=links[[var]](new_data[,var])
    sde_ci_df$X=links[[var]](new_data[,var])
    
    #labels for the plot
    labels=paste("A",1:n_id,sep="")
    
    p=ggplot()+geom_line(aes(X, par_estimates,col=ID),data=sde_par_df) +
      geom_ribbon(data=sde_ci_df,aes(x=X,ymin=lowpar,ymax=uppar,fill=ID),alpha=0.2)+
      scale_color_discrete(name="ID",labels=labels)+
      scale_fill_discrete(name="ID",labels=labels)+
      xlab(xlabels[[var]])+ylab(par)
    
    
    if (!(is.null(baseline))) {
      baseline_par=baseline$par(t = "all",new_data=new_data)
      baseline_CI=baseline$CI_pointwise(t="all",new_data=new_data)
      baseline_par_df=data.frame(ID=new_data$ID,Cov=new_data[,var],par_baseline = baseline_par[, par])
      baseline_par_df$X=links[[var]](new_data[,var])
      p=p+geom_line(aes(X, par_baseline,col=ID),data=baseline_par_df,linetype="dashed")
      }
    
    if (save) {
      
      if (!(dir.exists(model_name))) {
        dir.create(model_name)
      }
      
      ggsave(paste(paste("me",model_name,par,var,sep="_"),".png"),plot=p,width=10,height=5,path=model_name)
    }
    
    }
}




plot_fe_par_3D=function(baseline=NULL,sde,model_name,par,npoints=50,xmins=list(),xmaxs=list(),links=list(),
                        xlabels=list(),all_coeff_re_post,all_coeff_fe_post,level=0.95,save=TRUE) {
  
  #############################################################################################
  # 3D plot of a parameter when there are two non orthogonal covariates
  
  
  
  ###########################################################################################
  data=sde$data()
  
  #get the formula for the parameter of interest
  formulas=sde$formulas()
  form=formulas[[par]]
  
  #number of parameters
  n_par=length(names(formulas))
  
  #get all covariables in the model
  get_variables=function(formula) 
  {
    #list of strings for each term in the formula
    terms=colnames(attr(terms(formula),which="factors"))
    variables=unlist(lapply(terms, function(term) sub("s\\(([^,]+).*", "\\1", term)))
  }
  all_vars=unique(unlist(lapply(formulas,get_variables)))
  
  #get index of the parameter of interest in the vector of all parameters
  j=which(names(formulas)==par)[[1]]
  
  #list of strings for each term in the formula for the parameter of interest
  terms=colnames(attr(terms(form),which="factors")) 
  
  #fixed effects covariates/variables in the formula for the parameter of interest
  vars=get_variables(form) 
  vars=vars[vars!="ID"]
  
  #name of the two covariates
  var1=vars[1]
  var2=vars[2]
  print(var1)
  print(var2)
  
  #initialize data frame with covariate values 
  new_data<- data.frame(matrix(0, nrow = npoints*npoints, ncol = length(all_vars)))
  
  # Assign column names to the data frame
  colnames(new_data) <- all_vars
  
  
  #define the grid of covariates
  cov_values1=seq(from=xmins[[var1]],to=xmaxs[[var1]],length.out=npoints)
  cov_values2=seq(from=xmins[[var2]],to=xmaxs[[var2]],length.out=npoints)
  print(cov_values1)
  print(cov_values2)
  grid<- expand.grid(cov_values1,cov_values2)
  new_data[,vars]=grid
  
  #put any admissible value in the column ID
  new_data$ID=rep(data$ID[1],length(grid[,1]))
  
  
  #extract splines degrees of freedom
  str1=paste(par,paste("s(",var1,")",sep=""),sep=".")
  q1<- as.integer(sde$terms()$ncol_re[str1])+1
  
  str2=paste(par,paste("s(",var2,")",sep=""),sep=".")
  q2<- as.integer(sde$terms()$ncol_re[str2])+1
  
  #q1 <- as.integer(str_extract(terms[1], "(?<=k\\s=\\s)\\d+"))
  #q2 <-as.integer(str_extract(terms[2], "(?<=k\\s=\\s)\\d+"))
  
  #get model matrices 
  mats=sde$make_mat(new_data=new_data)
  X_fe=mats$X_fe
  X_re=mats$X_re
  
  #extract fixed effects coeff names
  fe_names1=paste(rep(paste(par,".s(",vars[1],").",sep=""),q1-1),1:(q1-1),sep="")
  fe_names2=paste(rep(paste(par,".s(",vars[2],").",sep=""),q1-1),1:(q1-1),sep="")
  fe_names=c(fe_names1,fe_names2)
  
  #keep coeffs and columns matching those names in the random effects
  coeff_re_names=rownames(sde$coeff_re())
  index=which(coeff_re_names %in% fe_names)
  X_re=X_re[,index]
  coeff_fe=sde$coeff_fe()
  coeff_re=sde$coeff_re()[index,]
  
  #linear predictor
  lp=X_fe%*%coeff_fe+X_re%*%coeff_re
  lp=matrix(lp,ncol=n_par)
  
  #parameters on response scale
  par_mat <- matrix(sapply(1:ncol(lp), function(i) {sde$invlink()[[i]](lp[,i])}), ncol = ncol(lp))
  colnames(par_mat) <- names(sde$invlink())
  
  #dataframe for fixed effects estimations
  est=as.data.frame(par_mat[,par])
  est$var1=new_data[,var1]
  est$var2=new_data[,var2]
  est$X1=links[[var1]](new_data[,var1])
  est$X2=links[[var2]](new_data[,var2])
  colnames(est)=c("par",var1,var2,"X1","X2")
  
  
  p <- plot_ly(est,type = "scatter3d",mode="markers",
               x = ~X1,y = ~X2,z=~par,color=~par)%>%
    layout(title=paste(par,"estimations"),scene=list(xaxis = list(title = xlabels[[var1]],showgrid = F),
    yaxis = list(title = xlabels[[var2]],showgrid = F),zaxis = list(title = par))) 
  
  #if we are given posterior draws of the parameters, we can plot the CIs
  if (!(is.null(all_coeff_fe_post)) & !(is.null(all_coeff_re_post))) {
  
    #confidence intervals for the population smooth
    coeff_re_post=all_coeff_re_post[,index]
    coeff_fe_post=all_coeff_fe_post
    
    n = nrow(X_fe)/n_par
    npost=length(all_coeff_fe_post[,1])
    
    #Get corresponding SDE parameters
    post_par <- array(sapply(1:npost, function(i) {
      sde$par(t = "all", 
                X_fe = X_fe, X_re = X_re, 
                coeff_fe =coeff_fe_post[i,],
                coeff_re = coeff_re_post[i,],
                resp = TRUE)  
    }), dim = c(n, n_par, npost))
    
    dimnames(post_par)[[2]] <- names(sde$invlink())
    
    alpha <- (1 - level)/2
    sde_CI <- apply(post_par, c(1, 2), quantile, probs = c(alpha, 1 - alpha))
    
    sde_CI <- aperm(sde_CI, c(3, 1, 2))
    dimnames(sde_CI)[[2]] <- c("low", "upp")
    
    # Data frame for CIs
    sde_ci_df <- data.frame(var1=new_data[,var1],var2=new_data[,var2],
                              lowpar = sde_CI[par, "low",],
                              uppar= sde_CI[par, "upp",],par=par_mat[,par])
    #add link variables
    sde_ci_df$X1=links[[var1]](new_data[,var1])
    sde_ci_df$X2=links[[var2]](new_data[,var2])
    
    # Create the initial 3D scatter plot
    
    
    p <- plot_ly(sde_ci_df,type = "scatter3d",mode="markers",
                 x = ~X1,y = ~X2,z=~par,color=~par)
    
      
    # Add the lower and higher estimation traces
    p %>% add_trace(x = ~X1,y = ~X2,z=~lowpar,
                    name = paste(par,"lower","estimation",sep=" "),
                    mode = "markers") %>%
        
    add_trace(x = ~X1,y = ~X2,z=~uppar,
                name = paste(par,"higher","estimation",sep=" "),
                mode = 'markers')
  
    
    # Layout adjustments for the 3D scene
    p %>%layout(title=paste(par,"estimations"),scene=list(xaxis = list(title = xlabels[[var1]],showgrid = F),       
          yaxis = list(title = xlabels[[var2]],showgrid = F),zaxis = list(title = par)))     
  }
  
  if (save) {
    
    if (!(dir.exists(model_name))) {
      dir.create(model_name)
    }
    
    # Save the plot as an HTML file
    file_path <- file.path(model_name, paste("fe", model_name, par, var1, var2, ".html", sep = "_"))
    saveWidget(p, file =file_path)
  }
  
  return (p)
}



####################### FUNCTION TO GET BEST INITIAL PARAMETERS ##########################


get_best_par0=function(sde,par0_df,fixpar=NULL) {
  
  formulas <- sde$formulas()
  data=sde$data()
  response=sde$response()
  type=sde$type()
  
  n=length(par0_df[,1])
  
  minAIC=Inf
  best_par0=c(0,0,0,0)
  
  for (i in 1:n) {
    par0 <- as.numeric(par0_df[i,])
    print(par0)
    baseline <- SDE$new(formulas = formulas,data = data,
                        type = type,response = response,par0 = par0,fixpar=fixpar)
    baseline$fit()
    
    AIC=baseline$AIC_marginal()
    #print(AIC)
    if (!is.nan(AIC) & AIC < minAIC) {
      best_par0=par0
      minAIC=AIC
    }
  }
  return (best_par0)
  
}

######################################################## FUNCTIONS TO CHECK MODELS ESTIMATIONS #############################################################################


check_fn=function(data) {
  
  res=c()
  
  for (id in unique(data$ID)) {
    
    #data for specific id
    data_id=data[data$ID==id,]
    n=length(data_id[,1])
    
    #time steps
    deltat=data_id[2:n,"time"]-data_id[1:(n-1),"time"]
    
    #empirical velocity 
    vx=(data_id[2:n,"x"]-data_id[1:(n-1),"x"])/deltat
    vy=(data_id[2:n,"y"]-data_id[1:(n-1),"y"])/deltat
    vnorm=sqrt(vx^2+vy^2)
    
    #mean and sd of velocity
    vmean=mean(vnorm)
    vsd=sqrt((sum((vnorm-vmean)^2))/(n-1))
    
    #angles between consecutive velocities
    thetas=c()
    for (i in 1:(n-2)) {
      v1=c(vx[i],vy[i])
      v2=c(vx[i+1],vy[i+1])
      thetas=c(thetas,signed_angle(v1,v2))
    }
    
    theta=mean(thetas)
  
    
    
    res=c(res,vmean,vsd,theta)
  }
  return (res)
  
}


plot_checks=function(check_stats,model_name) {
  
  n_id=6
  n_stats=3
  
  stats=check_stats$stats
  
  obs_stat=check_stats$obs_stat
  
  # Data frames for ggplot (observed and simulated values)
  df_obs <- as.data.frame.table(as.matrix(obs_stat))
  colnames(df_obs) <- c("stat", "sim", "val")
  df <- as.data.frame.table(stats)
  colnames(df) <- c("stat", "sim", "val")
  
  index=seq(from=1,to=length(df[,1]),by=n_stats)
  df_mean=df[index,]
  df_sd=df[index+1,]
  df_theta=df[index+2,]
  df_obs_mean=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats),]
  df_obs_sd=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats)+1,]
  df_obs_theta=df_obs[seq(from=1,to=n_id*n_stats,by=n_stats)+2,]
  
  # Create plot
  names=c(rep(paste("A",1:n_id,sep=""),each=n_stats))
  names(names)=paste("statistic",1:(n_id*n_stats),sep=" ")
  
  pmean <- ggplot(df_mean, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_mean,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  psd <- ggplot(df_sd, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_sd,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  ptheta <- ggplot(df_theta, aes(val)) + 
    geom_histogram(bins = 20, aes(y=..density..), col = "white", 
                   bg = "lightgrey", na.rm = TRUE) +
    facet_wrap(~stat, scales = "free",labeller=as_labeller(names)) + 
    geom_vline(aes(xintercept = val), data = df_obs_theta,col="red") +
    theme_light() +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black"))
  
  
  if (!(dir.exists(model_name))) {
    dir.create(model_name)
  }
  
  ggsave(paste(paste("mean_velocity_check",model_name,sep="_"),".png",sep=""),plot=pmean,path=model_name)
  ggsave(paste(paste("sd_velocity_check",model_name,sep="_"),".png",sep=""),plot=psd,path=model_name)
  ggsave(paste(paste("mean_angle_check",model_name,sep="_"),".png",sep=""),plot=ptheta,path=model_name)
}




check_fe_estimations=function(sde,model_name,links=list(),xmins=list(),xmaxs=list(),xlabels=list(),npost=1000,level=0.95) {
  
  #simulate a track without measurement errors from a fitted model
  #then fit the same model on this track 
  #and plot the estimated fixed effects with confidence interval vs the actual population smooth
  
  
  #get data
  data=sde$data()
  
  # get formulas and parameters
  formulas=sde$formulas()
  pars=names(formulas)
  n_par=length(formulas)
  
  #fixed coeffs (offsets)
  map=sde$map()
  
  #initial value of parameters 
  par0=sde$par()
  
  #simulate from the fitted model
  sde_sim=sde$simulate(data=data)
  
  #save plot of the sampled track
  plot_sample=ggplot(sde_sim)+ geom_point(mapping=aes(x, y,col=ID),size= 1,alpha=0.8) +
    scale_color_viridis_d(labels=c("A1","A2","A3","A4","A5","A6"))
  ggsave(paste(paste("check_fe",model_name,"sample",sep="_"),".png",sep=""),plot=plot_sample,width=10,height=5)
  
 
  #fit a new model on the simulated tracks
  sde_check=SDE$new(formulas = formulas,data = sde_sim,type = "sde",response = sde$response(),par0 = par0,map=map)
  sde_check$fit()
  
  #draw coeff from posterior distribution once and for all 
  # to get confidence intervals later
  postcoeff=sde_check$post_coeff(npost)
  all_coeff_re_post=postcoeff$coeff_re
  all_coeff_fe_post=postcoeff$coeff_fe
  
  #get all covariables in the model
  all_vars <- unique(unlist(lapply(formulas, function(formula) all.vars((formula)))))
  
  #Set link function to identity if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(links))) {
      link=(\(x) x)
      links[[var]]=link
    }
  }
  
  #Set xmins and xmaxs to min and max of observed covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xmins))) {
      xmin=min(data[,var])
      xmins[[var]]=xmin
    }
    if (!(var %in% names(xmaxs))) {
      xmax=max(data[,var])
      xmaxs[[var]]=xmax
    }
  }
  
  #Set xlabels to name of covariates if not given
  for (var in all_vars[all_vars!="ID"]) {
    if (!(var %in% names(xlabels))) {
      xlabels[[var]]=var
    }
  }
  
  
  #loop over the formulas of each parameter 
  for (j in seq_along(formulas)) {
    
    #formula for the parameter j
    form=formulas[[j]] 
    
    #name of parameter j
    par=pars[j] 
    print(par)
    
    #fixed effects covariates/variables in the formula
    vars=all.vars(form) 
    vars=vars[!(vars %in% c("ID","par0","n"))]
    
    #number of covariates
    ncovs=length(vars) 
      
    #if there is only one fixed effect covariable
    if (ncovs==1){
      
      link=links[[1]]
      xlabel=xlabels[[1]]
      xmax=xmaxs[[vars[1]]]
      xmin=xmins[[vars[1]]]
      
      
      p_check=plot_fe_par_2D(baseline=NULL,sde_check,model_name,par,npoints=200,xmin,xmax,link,xlabel,
                     all_coeff_re_post,all_coeff_fe_post,level,save=FALSE)
      p=plot_fe_par_2D(baseline=NULL,sde,model_name,par,npoints=200,xmin,xmax,link,xlabel,
                       all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level,save=FALSE)
      
      combined_plot <- ggplot() +
        geom_line(data = p_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
        geom_ribbon(data = p_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
        geom_line(data = p$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
        scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
        xlab(xlabels[[var]])+labs(color = "")
      
      if (!(dir.exists(model_name))) {
        dir.create(model_name)
      }
      
      ggsave(paste(paste("check_fe",model_name,par,var,sep="_"),".png",sep=""),plot=combined_plot,width=10,height=5,path=model_name)
      
      
    }
    
    # if there are two fixed effects covariables
    else if (ncovs==2){
      
      var1=vars[[1]]
      var2=vars[[2]]
      #if the covariates are orthogonal
      if (are_orthogonals(data,var1,var2)) {
        
        list_plots_check=plot_fe_par_2D_ortho(baseline=NULL,sde_check,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                             all_coeff_re_post,all_coeff_fe_post,level,save=FALSE)
        list_plots=plot_fe_par_2D_ortho(baseline=NULL,sde,model_name,par,npoints=200,xmins,xmaxs,links,xlabels,
                                        all_coeff_re_post=NULL,all_coeff_fe_post=NULL,level,save=FALSE)
        p1_check=list_plots_check[[1]]
        p2_check=list_plots_check[[2]]
        p1=list_plots[[1]]
        p2=list_plots[[2]]
        
        combined_plot1 <- ggplot() +
          geom_line(data = p1_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
          geom_ribbon(data = p1_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
          geom_line(data = p1$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
          scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
          xlab(xlabels[[var1]])+labs(color = "")
        
        combined_plot2 <- ggplot() +
          geom_line(data = p2_check$data, aes(x =X, y =par, color = "Fixed effects estimate"), size = 1) +
          geom_ribbon(data = p2_check$data,aes(x=X,ymin=lowpar,ymax=uppar),fill="red",alpha=0.2)+
          geom_line(data = p2$data, aes(x =X, y = par, color = "True fixed effects"), size = 1) +
          scale_color_manual(values = c("Fixed effects estimate" = "red", "True fixed effects" = "blue"))+
          xlab(xlabels[[var2]])+labs(color = "")
        
        if (!(dir.exists(model_name))) {
          dir.create(model_name)
        }
        
        ggsave(paste(paste("check_fe",model_name,par,var1,sep="_"),".png",sep=""),
               plot=combined_plot1,width=10,height=5,path=model_name)
        ggsave(paste(paste("check_fe",model_name,par,var2,sep="_"),".png",sep=""),
               plot=combined_plot2,width=10,height=5,path=model_name)
        
        
      }
      
      #else, plot in 3D
      else {
        
        p=plot_fe_par_3D(baseline,sde,model_name,par,npoints=200,xmins=xmins,xmaxs=xmaxs,links=links,
                       xlabels=xlabels,all_coeff_re_post,all_coeff_fe_post,level=0.95,save=FALSE)
      }
    }
  }
}




nearest_shore_point=function(point,coastline) {
  #find the nearest shorepoints on the coastline
  
  # Initialize variables to store the minimum distance and nearest point
  min_distance <- Inf
  nearest_point <- NULL
  
  
  # Iterate over each polygon in coastline_utm
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



simulate_constrained_ctcrw=function(z0,ctcrw,data,coastline) {
  
  ###############################################
  
  #z0 : initial point of the simulation
  #ctcrw : fitted ctcrw SDE
  #coastline : coastline data as UTM polygon
  
  #Simulate from a 2D CTCRW whose parameter mu depends
  #on the (normalized) normal vectors between the animal and the shore
  
  ###############################################"
  
  par <- ctcrw$par(new_data = data)
  
  # Loop over dimensions
  n_dim <- length(ctcrw$response())
  if(length(z0) < n_dim) {
    z0 <- rep(z0, n_dim)
  }
  # Initialize vector of simulated observations
  obs <- as.data.frame(matrix(rep(c(NA,NA), nrow(data)),nrow=nrow(data)))
  
  # Loop over IDs
  for(id in seq_along(unique(data$ID))) {
    
    # Get relevant rows of data
    print(id)
    ind <- which(data$ID == unique(data$ID)[id])
    dtimes <- diff(data$time[ind])
    sub_n <- length(ind)
    sub_par <- par[ind,]
    
    # hidden state
    alpha=matrix(c(z0[1],0,z0[2],0),nrow=4,ncol=1,byrow=TRUE)
    
    # Data has two columns for velocity and two for position
    sub_dat=data.frame("x1"=alpha[1],"v1"=alpha[2],"x2"=alpha[3],"v2"=alpha[4])
    
    # Unpack parameters
    tau <- sub_par[, n_dim + 1]
    nu <- sub_par[, n_dim + 2]
    beta <- 1/tau
    sigma <- 2 * nu / sqrt(tau * pi)
    
    # Loop over time steps
    mean <- rep(NA, 2)
    
    for(i in 2:sub_n) {
      
      #last location
      z=as.numeric(sub_dat[i-1,c("x1","x2")])
      #compute normal vector
      nearest_point=nearest_shore_point(st_point(1000*z),coastline)/1000
      nx=z[1]-nearest_point[1]
      ny=z[2]-nearest_point[2]
      Dshore=sqrt(nx^2+ny^2)
      normalized_nx=nx/Dshore
      normalized_ny=ny/Dshore
      
      #adjust covariate data
      new_data=data[i-1,]
      new_data$normalized_nx=normalized_nx
      new_data$normalized_ny=normalized_ny
      
      #get value of mu
      mu=ctcrw$par(new_data=new_data)[1,c("mu1","mu2")]
      mu=matrix(mu,ncol=1,nrow=2,byrow=TRUE)
      
      
      #covariance matrix
      var_xi=sigma[i-1]^2/beta[i-1]^2*(dtimes[i-1]+2/beta[i-1]*(exp(-beta[i-1]*dtimes[i-1])-1)+
                                         1/2/beta[i-1]*(1-exp(-2*beta[i-1]*dtimes[i-1])))
      var_zeta=sigma[i-1]^2/2/beta[i-1]*(1-exp(-2*beta[i-1]*dtimes[i-1]))
      cov=sigma[i-1]^2/2/beta[i-1]^2*(1+exp(-2*beta[i-1]*dtimes[i-1])-2*exp(-beta[i-1]*dtimes[i-1]))
      Qi=matrix(c(var_xi,cov,cov,var_zeta),nrow=2,byrow=TRUE)
      Qi_tilde=bdiag(replicate(2,Qi,simplify=FALSE))
      
      #drift part
      Bi=matrix(c(dtimes[i-1],0),nrow=2,ncol=1,byrow=TRUE)
      Bi_tilde=bdiag(replicate(2,Bi,simplify=FALSE))
      
      #random part
      eta=mvrnorm(1,mu =c(0,0,0,0), Sigma =Qi_tilde)
      
      #hidden states link matrix
      Ti=matrix(c(1,1/beta[i-1]*(1-exp(-beta[i-1]*dtimes[i-1])),0,exp(-beta[i-1]*dtimes[i-1])),nrow=2,byrow=TRUE)
      Ti_tilde=bdiag(replicate(2,Ti,simplify=FALSE))
      
      #state obs link matrix
      Z=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE)
      Z_tilde=bdiag(replicate(2,Z,simplify=FALSE))
      
      
      alpha=Ti_tilde%*%alpha+Bi_tilde%*%mu+eta
      sub_dat=rbind(sub_dat,c(alpha[1,1],alpha[2,1]+mu[1],alpha[3,1],alpha[4,1]+mu[2]))
      
    }
    # Only return location (and not velocity)
    sub_obs <- sub_dat[,c("x1","x2")]
    # Update observation vector
    obs[ind,] <- sub_obs
  }

  # Add simulated variable to data frame
  data[,ctcrw$response()] <- obs
  return (data)
}


#' @description Compute signed angle in [-pi,pi] that rotates first vector into second vector
#' 
#' @param point the location we want to tell if it is on land or not
#' @param land sf object (list of polygons) defining the land
#' @return boolean TRUE if it is on land, FALSE otherwise
#' 
signed_angle <- function(u, v) {
  # u, v: matrices of bivariate vectors with 2 columns
  u <- matrix(u, ncol = 2)
  v <- matrix(v, ncol = 2)
  if (nrow(u) != nrow(v)) stop("u and v must have the same number of 
                                  rows")
  result <- as.numeric(atan2(v[,2], v[,1]) - atan2(u[,2], u[,1]))
  ind1 <- which(result > pi)
  ind2 <- which(result <= -pi)
  result[ind1] <- result[ind1] - 2*pi
  result[ind2] <- result[ind2] + 2*pi
  return(result) 
} 
