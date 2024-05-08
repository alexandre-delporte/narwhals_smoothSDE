
#####################              CHECK GOODNESS OF FIT OF BASELINE MODEL            ################################



check_baseline1=baseline1$check_post(check_fn, n_sims =200, silent = TRUE)
plot_checks(check_baseline1,"baseline1")

check_baseline2=baseline2$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline2,"baseline2")

check_baseline3=baseline3$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline3,"baseline3")

check_baseline4=baseline4$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline4,"baseline4")

check_baseline5=baseline5$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline5,"baseline5")

check_baseline6=baseline6$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline6,"baseline6")

check_baseline7=baseline7$check_post(check_fn, n_sims =500, silent = TRUE)
plot_checks(check_baseline7,"baseline7")




