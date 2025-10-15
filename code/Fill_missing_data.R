fill_missing_data <- function(path_df)
{
     this_t <- path_df$POSITION_T
     this_x <- path_df$POSITION_X
     this_y <- path_df$POSITION_Y
  this_t = round(this_t,digits=4)
  this_dt = round(min(diff(this_t)),digits=4)
  diff_t = round(diff(this_t),digits=4)
  this_new_t = seq(this_t[1],this_t[length(this_t)],by=this_dt)
  this_new_x = 0*this_new_t
  this_new_y = 0*this_new_t
  
  sigma_est_x = sd(diff(this_x))/2
  sigma_est_y = sd(diff(this_y))/2
  
  find_next_good = function(x,idx) {
    i=idx
    while(x[i] == 0){
      i = i+1
    }
    return(i)
  }
  
  for (tt in 1:length(this_t)) {
    this_new_idx = which(abs((this_new_t - this_t[tt]))<0.0005)
    this_new_x[this_new_idx] = this_x[tt]
    this_new_y[this_new_idx] = this_y[tt]
  }
  
  fill_idx_x = which(this_new_x == 0)
  fill_idx_y = which(this_new_y ==0)
  
  last_good_idx = 1
  for (tt in 2:(length(this_new_t)-1)){
    #print(tt)
    if (this_new_x[tt] == 0 || this_new_y[tt]==0) {
      next_good_idx_x = find_next_good(this_new_x,tt)
      next_good_idx_y = find_next_good(this_new_y,tt)
      
      this_new_x[tt] = (this_new_t[next_good_idx_x] - this_new_t[tt])/
        (this_new_t[next_good_idx_x]-this_new_t[last_good_idx])*
        this_new_x[last_good_idx] +
        (this_new_t[tt] - this_new_t[last_good_idx])/
        (this_new_t[next_good_idx_x]-this_new_t[last_good_idx])*
        this_new_x[next_good_idx_x] + sigma_est_x*rnorm(1,0,1)
      
      this_new_y[tt] = (this_new_t[next_good_idx_y] - this_new_t[tt])/
        (this_new_t[next_good_idx_y]-this_new_t[last_good_idx])*
        this_new_y[last_good_idx] +
        (this_new_t[tt] - this_new_t[last_good_idx])/
        (this_new_t[next_good_idx_y]-this_new_t[last_good_idx])*
        this_new_y[next_good_idx_y] + sigma_est_y*rnorm(1,0,1)
    } else {
      last_good_idx = tt
    }
  }
  
  
  
  
  
  #points(this_new_t[replace_idx],this_new_x[replace_idx],col="red",pch=4)
  
  replace_idx_x = which(diff(this_new_x) > 5*sigma_est_x) + 1
  replace_idx_y = which(diff(this_new_y) > 5*sigma_est_y) + 1
  
  for (tt in replace_idx_x) {
    if((tt > 1&& tt<length(this_new_t)) || tt < length(this_new_idx)) {
      last_good_idx_x = tt-1
      next_good_idx_x = tt+1
      this_new_x[tt] = (this_new_t[next_good_idx_x] - this_new_t[tt])/
        (this_new_t[next_good_idx_x]-this_new_t[last_good_idx_x])*
        this_new_x[last_good_idx_x] +
        (this_new_t[tt] - this_new_t[last_good_idx_x])/
        (this_new_t[next_good_idx_x]-this_new_t[last_good_idx_x])*
        this_new_x[next_good_idx_x] + sigma_est_x*rnorm(1,0,1)
    }
  }
  
  for (tt in replace_idx_y) {
    if ((tt > 1&& tt<length(this_new_t))|| tt < length(this_new_idx)) {
      last_good_idx_y = tt-1
      next_good_idx_y = tt+1
      this_new_y[tt] = (this_new_t[next_good_idx_y] - this_new_t[tt])/
        (this_new_t[next_good_idx_y]-this_new_t[last_good_idx_y])*
        this_new_y[last_good_idx_y] +
        (this_new_t[tt] - this_new_t[last_good_idx_y])/
        (this_new_t[next_good_idx_y]-this_new_t[last_good_idx_y])*
        this_new_y[next_good_idx_y] + sigma_est_y*rnorm(1,0,1)
    }
  }
  plot(this_new_t,this_new_x,pch=20)
  points(this_new_t[fill_idx_x],this_new_x[fill_idx_x],col="red",pch=20)
  points(this_new_t[replace_idx_x],this_new_x[replace_idx_x],col="blue",pch=4)
  
  plot(this_new_t,this_new_y,pch=20)
  points(this_new_t[fill_idx_y],this_new_y[fill_idx_y],col="red",pch=20)
  points(this_new_t[replace_idx_y],this_new_y[replace_idx_y],col="blue",pch=4)
  
  return(list(t=this_new_t,
              x=this_new_x, 
              y=this_new_y))
}
