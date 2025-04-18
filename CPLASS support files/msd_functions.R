pathwise_MSD = function(path,use_pct = 0.2){
  x = path$path$x
  y = path$path$y
  t = path$path$t
  
  n = length(t)
  n_lag = floor(use_pct*n)
  t_lag = t[1:n_lag]
  msd = 0*t_lag
  for (tt in 1:n_lag){
    this_diff = rep(NA,n-tt)
    for (ttt in 1:(n-tt)){
      this_diff[ttt] = (x[ttt+tt] - x[ttt])^2 + (y[ttt+tt] - y[ttt])^2
    }
    msd[tt] = mean(this_diff)/4
  }
  output = tibble(
    t_lag = t_lag,
    msd = msd
  )
  return(output)
}

ensemble_MSD = function(msd_df) {
  num_paths = length((msd_df %>% distinct(id))$id)
  num_lags = length((msd_df %>% distinct(t))$t)
  
  e_msd = tibble()
  for (tt in (msd_df %>% distinct(t))$t){
    
    e_msd = bind_rows(e_msd,tibble(
      t = tt,
      e_msd = mean((msd_df %>% filter(t == tt))$msd)
    ))
  }
  return(e_msd)
}

plot_msd = function(msd_list,is_active,group_name = "") {
  emsd = ensemble_MSD(msd_list)
  full_msd_t = emsd$t_lag
  dt = full_msd_t[2]-full_msd_t[1]
  active_list = which(is_active ==  TRUE)
  not_active = which(is_active ==  FALSE)
  plot(1,type="l",xlim = c(dt,max(full_msd_t)), ylim = c(0.001,1),
       log='xy', xlab = "Time (s)", ylab = "Microns^2",
       main = paste0("Group MSD: ",group_name))
  msd_not_active = list()
  msd_active = list()
  count = 1
  for (i in not_active){
    lines(msd_list[[i]]$t_lag,msd_list[[i]]$msd,type="l",
          col="pink")
    msd_not_active[[count]] = msd_list[[i]]
    count=count+1
  }
  count = 1
  for (i in active_list){
    lines(msd_list[[i]]$t_lag,msd_list[[i]]$msd,type="l",
          col="darkgray")
    msd_active[[count]] = msd_list[[i]]
    count=count+1
  }
  emsd_not_active = ensemble_MSD(msd_not_active)
  emsd_active = ensemble_MSD(msd_active)
  lines(emsd_not_active$t_lag,emsd_not_active$e_msd,lwd=2,col="red")
  lines(emsd_active$t_lag,emsd_active$e_msd,lwd=2,col="black")
  lines(emsd$t_lag,emsd$e_msd,lwd=2,col="blue")
  legend("topleft",
         legend = c("Ensemble MSD",
                    "Pathwise MSD, Not Active",
                    "Ensemble MSD, Not Active",
                    "Pathwise MSD, Active",
                    "Ensemble MSD, Active"),
         col = c("blue","pink","red","gray","black"),
         lwd = c(2,1,2,1,2)
  )
  return(emsd)
}

