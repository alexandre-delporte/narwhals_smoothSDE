


#################################    POSTERIOR PREDICTIVE CHECKS  #########################################


check_response2=response2$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2,"response2")

check_response3=response3$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3,"response3")

check_response4=response4$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4,"response4")

check_response5=response5$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5,"response5")

check_response6=response6$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6,"response6")

check_response7=response7$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response7,"response7")

check_response8=response8$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response8,"response8")

check_response9=response9$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response9,"response9")

check_response10=response10$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response10,"response10")

check_response11=response11$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response11,"response11")

check_response12=response12$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response12,"response12")

check_response13=response13$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response13,"response13")





check_response2offset=response2offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response2offset,"response2offset")

check_response3offset=response3offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response3offset,"response3offset")

check_response4offset=response4offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response4offset,"response4offset")

check_response5offset=response5offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response5offset,"response5offset")

check_response6offset=response6offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response6,"response6offset")

check_response7offset=response7offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response7offset,"response7offset")

check_response8offset=response8offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response8offset,"response8offset")

check_response9offset=response9offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response9offset,"response9offset")

check_response10offset=response10offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response10offset,"response10offset")

check_response11offset=response11offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response11offset,"response11offset")

check_response12offset=response12offset$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response12offset,"response12offset")

check_response13offset=response13offste$check_post(check_fn, n_sims = 500, silent = TRUE)
plot_checks(check_response13offset,"response13offset")







#######################   SIMULATE DATA WITH FITTED PARAMETERS MODEL AND  TRY TO RECOVER IT #########################

# models without offset
check_fe_estimations(ctcrw=response2,model_name="response2",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3,model_name="response3",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4,model_name="response4",links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)





#models with offset
check_fe_estimations(ctcrw=response2offset,model_name="response2offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response3offset,model_name="response3offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response4offset,model_name="response4offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)


check_fe_estimations(ctcrw=response8offset,model_name="response8offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)

check_fe_estimations(ctcrw=response9offset,model_name="response9offset",
                     links=list("ExpShip"=(\(x) 1/x)),
                     xmins=list("ExpShip"=0.02),xmaxs=list("ExpShip"=0.2),
                     xlabels=list("ExpShip"="Distance to ship"),npost=1000,level=0.95)




