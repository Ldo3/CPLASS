#Running CPLASS for a single path
library(here)
library("Rlab")
library("matrixcalc")
library("pracma")
library("matrixcalc")
library("limSolve")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("foreach")
library("doParallel")
library(igraph)
library(readr)

CPLASS<- function(t,x, y, lambda_r=1/30, 
                  iter_max = 5000, burn_in=500, s_cap=1, 
                  gamma=1.01, speed_pen=TRUE, Diagnostic = FALSE, 
                  sd = NA, pen = "ssic")
{
     dt = min(diff(t))
     n <- length(t)

     if (speed_pen == TRUE) {
          speed_control <- 1
     } else {
          speed_control <- 0
     }

     # Auto-retry loop
     max_retries <- 20  # set a reasonable retry limit
     attempt <- 1
     pl = list()
     repeat {
          MCMC <- try(MHsearch(t, x, y, dt, lambda_r = lambda_r,
                            iter_max = iter_max, burn_in = burn_in,
                            s_cap = s_cap, gamma = gamma,
                            speed_control = speed_control, sd = NA, pen ="ssic"),
                      silent = TRUE)

          if (!inherits(MCMC, "try-error")) {
               # Success
               Final_MCMC = FinalMH(MCMC$cps_list)
               index.max = which.max(Final_MCMC$criterion_score)
               pl = piecewise_linear_con(t, x, y,
                                          Final_MCMC$this_cp[[index.max]],
                                          old_version = FALSE)
               if(Diagnostic == TRUE)
               {
                 info_table = MCMC$info_table
                 update_info = MCMC$update_info
               }
               break
          } else {
               message(paste("MHsearch failed, retrying... attempt", attempt))
               attempt = attempt + 1
               if (attempt > max_retries) {
                    message("MHsearch failed too many times. Returning pl = NULL.")
                    break
               }
          }
     }

     if(Diagnostic == TRUE){return(list(pl = pl,
                                        info_table = info_table,
                                        update_info = update_info))}else{return(pl)}

}

# MHsearch() : Metropolis-Hastings algorithm with tailored proposal functions for stochastics search
MHsearch<- function(t,x, y,dt, lambda_r=1/30, 
                    iter_max = 5000, burn_in=500, s_cap=1, 
                    gamma=1.01, speed_control=0, sd = NA, pen ="ssic")
{
     #N: sample size
     #r0: is the initial vector of switching process, this (N-2) dim vector contains 0s, 1s
     #dt: e.g. 0.05 for 20Hz, 0.1 for 10Hz, 1 for 100Hz, 0.02 for 50Hz
     #lambda_r: the rate in Bernoulli process used to propose a random vector of changepoints r
     #iter_max: number of interations in running Markov Chain Monte Carlo
     #burn_in: number of burn in steps
     #s_cap: a minimum speed that has no penalty
     #gamma: the power in the strengthed Schwarz Criterion form
     #speed_control: choose 0 for deactivate the speed penalty, choose 1 for activate the speed penalty


     # Initial values
     N = length(t)
     penalty_coef = max(4,log(N))^gamma # strengthened Schwarz Criterion
     r0 = q_new(lambda_r, dt = dt, N) #not include the starting and the ending points of the path
     indexcp= which(r0==1)
     cps0 = indexcp+1 #initial index of changepoint times
     initial_input = piecewise_linear_con(t, x, y, cps0)
     seg_vel0 = initial_input$path_segvel #initial segment speeds
     seg_time0= initial_input$path_segtime #initial segment durations
     u0 = initial_input$u_x #initial segment velocities w.r.t x-axis
     v0 = initial_input$v_y #initial segment velocities w.r.t y-axis

     RSS0 = initial_input$path_RSS #initial residual sum of square


     cps_list = list() #store the list of changepoints in all iterations
     FLAG = TRUE
     info_table = tibble(llh_cur = c(),
                         llh_new = c(),
                         score_cur = c(),
                         score_new = c(),
                         ncp_cur = c(),
                         ncp_new = c(),
                         pen_cur = c(),
                         pen_new = c(),
                         est_sigma_cur = c(),
                         est_sigma_new = c(),
                         ind_iteration = c(),
                         decision = c())

     # Metropolis-Hastings algorithm to update r
     pb <- txtProgressBar(min = 0, max = iter_max, style = 3)  # style=3 gives percentage
     for (count in 1:iter_max)
     {
          # print(count)
          # if((count %% 1000)==0){print(count)}
          k = length(v0) # number of segments

          ######## Update r | data <- MH ######

          ##### RUNNING METROPOLIS-HASTING ######
          u = runif(1)
          pp = proposal_function(u, r0, N, lambda_r, dt = dt)
          r_prop = pp$r_prop
          status = pp$status

          # Log Acceptance function CS<-function(t,x,y,cps)

          A = 1
          l_prop = CS(t, x, y, r_prop, penalty_coef, 
                      s_cap, speed_control, sd = sd, 
                      pen = pen, gamma = gamma)
          l_cur = CS(t, x, y, r0, penalty_coef, s_cap, 
                     speed_control, sd = sd, pen = pen, 
                     gamma = gamma)
          if (l_prop$logical ==FALSE || l_cur$logical == FALSE){
               A =0
          }else if(l_prop$s ==0){A=0}else{
               logA = l_prop$s+pproposal(u, r0, N, lambda_r, r_prop, 1-status, dt = dt)-l_cur$s-pproposal(u, r_prop, N, lambda_r, r0, status, dt = dt)
          }


          if((logA >=0  || runif(1) < exp(logA)) & A!=0 ){r1 = r_prop}else{r1 = r0
          }

          index1 = which(r1 ==1)
          cps1 = index1 + 1 #proposed changepoints
          np = piecewise_linear_con(t,x,y,cps1) #fit piecewise linear to this proposed changes
          seg_vel1 = np$path_segvel #the inferred segment speeds
          seg_time1= np$path_segtime #the inferred segment durations
          u1 = np$u_x # the inferred segment velocities x-axis
          v1 = np$v_y # the inferred segment velocities y-axis
          RSS1 = np$path_RSS
          criterion_score = CS(t ,x ,y ,r1 ,penalty_coef ,s_cap, speed_control, sd =sd, pen = pen, gamma =gamma)$s


          cps_list[[count]] = list(r = r1,
                                   u = u1,
                                   v = v1,
                                   this_cp = cps1,
                                   seg_vel= np$path_segvel,
                                   seg_time = np$path_segtime,
                                   eta = np$path_segeta,
                                   criterion_score = criterion_score,
                                   RSS = RSS1,
                                   r_prop = r_prop
          )
          if(any(r1!=r0)){decision = 1}else{decision = 0}
          info_t = tibble(llh_cur = l_cur$llh,
                              llh_new = l_prop$llh,
                              score_cur = l_cur$s,
                              score_new = l_prop$s,
                              ncp_cur = sum(r0),
                              ncp_new = sum(r_prop),
                              logA = logA,
                              pen_cur = l_cur$p,
                              pen_new = l_prop$p,
                              est_sigma_cur = sqrt(RSS0/(2*length(t))),
                              est_sigma_new = sqrt(RSS1/(2*length(t))),
                              ind_iteration = count,
                              decision = decision)
          info_table = bind_rows(info_table, info_t)

          #Record updating r0, c0, u0, v0
          r0 = r1
          u0 = u1
          v0 = v1
          seg_vel0 = seg_vel1
          seg_time0 = seg_time1
          RSS0 = RSS1

          setTxtProgressBar(pb, count)
     }
     close(pb)
     # print(paste("Reject rate is = ",reject/iter_max, sep="" ))
     # print(paste("Reject qd rate is = ",count_reject_qd/count_qd, sep="" ))
     # print(paste("Reject q_new rate is =", count_reject_q_new/count_q_new, sep = ""))
     # print(paste("Reject q_shift rate is = ", count_reject_q_shift/count_q_shift,sep = ""))
     #
     info_table <- info_table %>% mutate_all(~replace(., is.na(.), 0))
     it = info_table %>% slice(-(1:burn_in))
     update_info = tibble(accept_rate = sum(it$decision)/ length(it$decision),
                          reject_rate = 1-sum(it$decision)/ length(it$decision))

     return(list(cps_list = cps_list[-(1:burn_in)],
                 info_table = it,
                 update_info = update_info))
}



# Proposal functions:

#1. Independent switch point process
log_q_new_pmf <- function(lambda_r, dt, N, r)
{
     K_r = sum(r) #number of change points
     p = K_r * log(1-exp(-lambda_r*dt)) - (N-2-K_r)*(lambda_r*dt)
     return(p) #output is a number
}

q_new <- function(lambda_r, dt, N)
{
     p = 1-exp(-lambda_r*dt)
     return(rbern(N-2,p))
}


#2. The creation or extinction of a change point
log_q_bd_pmf <- function(r,status)
{
     N_r= length(r) #N_r=N-2 (= n-1 on the paper)
     K_r = sum(r) #number of change points
     if(status == 0){ p = 1/(2*K_r)}

     if(!any(r==1) && status == 0){p=0}

     if(status == 1){p= 1/(2*(N_r-K_r))}

     if(status == 1 && !any(r==0)){p= 0}
     return(log(p)) #output is a number
}


q_bd <- function(r_cur)
{
     r_prop = r_cur
     N_r = length(r_cur)
     s = sample(c(1:N_r), size = 1)
     if (!any(r_cur==0)||!any(r_cur==1)){
          for(i in 1:length(r_cur))
          {
               if(i==s){r_prop[i]=1-r_cur[i]
               status = r_prop[i]
               }else{r_prop[i]=r_cur[i]}
          }
          }else{
          status =sample(c(0,1),size=1) # 1 for adding a changepoint, 0 for removing a changepoint
          if(status==1){s=sample(which(r_cur==0),size=1)
          r_prop[s]=1
          }else{
               s= sample(which(r_cur==1),size=1)
               r_prop[s]=0}
          }

     return(list(r_prop=r_prop,
                 status = status)) #return a proposal vector r
}
#3. a location shift of a single change point
q_shift <- function(r_cur)
{
     ones_idx <- which(r_cur == 1)
     zeros_idx <- which(r_cur == 0)

     if (length(ones_idx) == 0 || length(zeros_idx) == 0) {
          stop("Error: Need at least one 1 and one 0 in the vector.")
     }

     # print(paste("ones at:", paste(ones_idx, collapse = ",")))
     # print(paste("zeros at:", paste(zeros_idx[1:5], collapse = ","), "..."))

     if(length(ones_idx) ==1){s = ones_idx}else{s <- sample(ones_idx, size = 1)}
     if(length(zeros_idx)==1){ss = zeros_idx}else{ss <- sample(zeros_idx, size = 1)}


     stopifnot(r_cur[s] == 1)
     stopifnot(r_cur[ss] == 0)

     r_prop <- r_cur
     r_prop[c(s, ss)] <- 1 - r_cur[c(s, ss)]

     # cat("Selected s (1→0):", s, "value:", r_cur[s], "\n")
     # cat("Selected ss (0→1):", ss, "value:", r_cur[ss], "\n")
     # cat("Before:", sum(r_cur), " ones\n")
     # cat("After :", sum(r_prop), " ones\n")

     return(r_prop)
}



# #4. The creation or extinction of a segment
q_as_pmf <- function(r)
{
     N_r= length(r) #N_r=N-2 (= n-1 on the paper)
     K_r = sum(r) #number of change points
     #idx_add=sample(c(1:(K_r+1)),1)
     cps = which(r==1)+1
     cp = append(cps,1)
     cp = append(cp,N_r+1)
     cp = sort(cp)
     d=diff(cp)
     com=c()
     for(i in 1:length(d)){
          if(d[i]<2){com[i]=0}else{com[i]=(d[i]-1)*(d[i]-2)}
     }
     p= 1/(2*((N_r-K_r)*(N_r-K_r-1)))*sum(com)
     return(p) #output is a number
}

q_as <- function(r_cur)
{
     r_prop=r_cur
     N_r= length(r_cur) #N_r=N-2 (= n-1 on the paper)
     K_r = sum(r_cur) #number of change points
     idx_add=sample(c(1:(K_r+1)),1)
     cp_cur = which(r_cur==1)+1
     cp = append(cp_cur,1)
     cp = append(cp,N_r+2)
     cp = sort(cp)
     if(cp[idx_add+1]-cp[idx_add]>2){new = sample(c(cp[idx_add]:cp[idx_add+1]),2)
     cp_prop=unique(sort(append(cp_cur,new)))
     r_prop[cp_prop-1]=1
     }
     return(r_prop) #return a proposal vector r
}

# Delete 2 consecutive cps
q_ds_pmf <- function(r)
{
     N_r= length(r) #N_r=N-2 (= n-1 on the paper)
     K_r = sum(r) #number of change points
     if (K_r >= 2){p = 1/K_r}else{
          p = 0
     }
     return(p) #output is a number
}

q_ds <- function(r_cur)
{
     r_prop=r_cur
     if(r_cur[length(r_cur)]==1&&length(which(r_cur==1))<=2){
          s=which(r_cur==1)[-length(which(r_cur==1))]}else{
               s= sample(which(r_cur==1)[-length(which(r_cur==1))],size=1)
          }
     r_prop[s]=0
     r_prop[which(r_cur==1)[which(which(r_cur==1)==s)+1]]=0
     return(r_prop) #return a proposal vector r
}


# Delete and add a segment
q_bd2<-function(r_cur)
{
     r_prop=r_cur
     if(!any(r_cur==1)){
          s=sample(c(1:length(r_cur)),2)
          r_prop[s]=1
          status = 1
     }else if(!any(r_cur==0)){
          s=sample(c(1:length(r_cur)),1)
          r_prop[s]=0

          if(s==length(r_cur)){r_prop[s-1]=0}else{
               r_prop[s+1]=0}

          status = 0
     }else if(length(which(r_cur==1))==1){
          r_prop=q_as(r_cur)
          status = 1
     }else{
          status = sample(c(0,1),1)
          if(status==1){#adding a segment
               r_prop=q_as(r_cur)}else{r_prop=q_ds(r_cur)}
          }

     return(list(r_prop = r_prop,
                 status = status))
}

log_q_bd2_pmf<-function(r, status)
{
     if(status == 0){p = q_ds_pmf(r)}else{p = q_as_pmf(r)}

     return(log(p))
}
# Final proposal function
proposal_function<-function(u, r_cur, N, lambda_r, dt)
{
     if (!any(r_cur==0)||!any(r_cur==1)){
          if (0<=u && u<=1/4){
               r_prop = q_new(lambda_r, dt, N)
               status = 2
               }else if(1/4<u &&u<=1/2){
                    pp = q_bd(r_cur)
                    r_prop = pp$r_prop
                    status = pp$status
                    }else{
                         pp = q_bd2(r_cur)
                         r_prop = pp$r_prop
                         status = pp$status
                         }
          }else{
               if (0<=u && u<=1/4){
                    r_prop = q_new(lambda_r, dt, N)
                    status = 2
                    }else if(1/4<u &&u<=3/8){
                         pp = q_bd(r_cur)
                         r_prop = pp$r_prop
                         status = pp$status
                         }else if(3/8<u && u<=1/2){
                              pp = q_bd2(r_cur)
                              r_prop = pp$r_prop
                              status = pp$status
                              }else{
                                   r_prop = q_shift(r_cur)
                                   status = 3
                                   }
          }
     return(list(r_prop = r_prop,
                 status = status))
}

### PDF of the proprosal function
pproposal<-function(u, r, N, lambda_r, r_given, status, dt)
{
     if (!any(r_given==0)||!any(r_given==1)){
          if (0<=u && u<=1/4){
               logp= log_q_new_pmf(lambda_r,dt, N, r)
               }else if(1/4<u &&u<=1/2){
                    logp = log_q_bd_pmf(r_given,status)
                    }else{
                         logp = log_q_bd2_pmf(r_given,status)
                         }
          }else{
               if (0<=u && u<=1/4){
                    logp = log_q_new_pmf(lambda_r, dt, N, r)
                    }else if(1/4<u &&u<=3/8){
                         logp = log_q_bd_pmf(r_given,status)
                         }else if(3/8<u && u<=1/2){
                              logp = log_q_bd2_pmf(r_given,status)
                              }else{
                                   logp = log(1/(N-2))
                              }
          }
     return(logp)
}



# Log-likelihood Penalty: See PENALTY.R file 



# The following function is used for reading results after running MH algorithm
FinalMH <- function(bcp_new)
{
     r=list()
     n = length(bcp_new)
     for (i in 1:n)
     {
          r[[i]]=bcp_new[[i]]$r
     }

     u=list()
     for (i in 1:n)
     {
          u[[i]]=bcp_new[[i]]$u
     }

     v=list()
     for (i in 1:n)
     {
          v[[i]] =bcp_new[[i]]$v
     }


     seg_vel=list()
     for (i in 1:n)
     {
          seg_vel[[i]] =bcp_new[[i]]$seg_vel
     }

     seg_time=list()
     for (i in 1:n)
     {
          seg_time[[i]] =bcp_new[[i]]$seg_time
     }

     this_cp =list()
     for (i in 1:n)
     {
          this_cp[[i]] = bcp_new[[i]]$this_cp
     }

     eta = c()
     for (i in 1:n)
     {
          eta[i] = bcp_new[[i]]$eta
     }

     loglikelihood = c()
     for (i in 1:n)
     {
          loglikelihood[i] =bcp_new[[i]]$loglikelihood
     }

     RSS = c()
     for (i in 1:n)
     {
          RSS[i] =bcp_new[[i]]$RSS
     }

     criterion_score = c()
     for (i in 1:n)
     {
          criterion_score[i] =bcp_new[[i]]$criterion_score
     }

     r_prop=list()
     n = length(bcp_new)
     for (i in 1:n)
     {
          r_prop[[i]]=bcp_new[[i]]$r_prop
     }
     num_cp =c()
     for (i in 1:n)
     {
          num_cp[i] = length(bcp_new[[i]]$this_cp)
     }
     return(list(r=r,
                 u = u,
                 v = v,
                 this_cp = this_cp,
                 seg_vel= seg_vel,
                 seg_time = seg_time,
                 eta = eta,
                 criterion_score = criterion_score,
                 RSS=  RSS,
                 num_cp= num_cp,
                 r_prop = r_prop
     ))
}



## PLA


# We will construct a code to due with piecewise linear regression needed
#to be continuous on the entire domain

safe_qrsolve <- function(X, y) {
  tryCatch({
    coef <- qr.solve(X, y)
    list(value = coef, logical = TRUE)
  }, error = function(e) {
    list(value = NA, logical = FALSE)
  })
}

#1. the first version is a simplify version to run inside the MH search
piecewise_linear_con<-function(t,x,y,this_cp, old_version=TRUE)
{
     #this_cp is the vector containing index of change points
     n = length(x) #number of data points
     m = length(this_cp) +1  # number of segment, m-1 is number of change points {1,...,m-1}
     # so we have m segments [0,c1],(c1,c2],...,(c_{m-1},n]
     this_cp = sort(this_cp)
     this_segvel = c()
     this_segtime = c()
     this_segtheta = c()  # Note this takes values between -pi/2 and pi/2
     this_segeta = c()
     vx_vec = c()
     vy_vec = c()

     if(length(this_cp)>0){
          # Contructing the regression matrix A

          A = matrix(1, n , m+1) #require (n>(m+1))
          # Update for the second column
          for (i in 1:n)
          {
               A[i,2] = t[i]
          }
          # Update from the third column to the (m+1)th column
          for (i in 1:n)
          {
               for (j in 3:(m+1))
               {
                    A[i,j] = (t[i]-t[this_cp[j-2]])*d_fn(i,this_cp[j-2])
               }
          }
     }else{
          A = matrix(1,n,2)
          for (i in 1:n)
          {
               A[i,2] = t[i]
          }
     }
     alphaS = safe_qrsolve(A, x)
     betaS= safe_qrsolve(A, y)

     #estimate for vector alpha and beta, sigma
     if ((alphaS$logical == FALSE || betaS$logical == FALSE)){return(list(logical=FALSE))}else{
          # alpha = inv(t(A) %*% A) %*% t(A) %*% x
          # beta = inv(t(A) %*% A) %*% t(A) %*% y
          alpha = alphaS$value
          beta = betaS$value

          alpha_intercept = alpha[1]
          beta_intercept = beta[1]


          x_plinear = A %*% alpha
          x_plinear = as.vector(x_plinear)
          y_plinear = A %*% beta
          y_plinear = as.vector(y_plinear)
          RSS =  t(x - A %*% alpha) %*% (x - A %*% alpha) + t(y - A %*% beta) %*% (y - A %*% beta)
          RSS = as.numeric(RSS)
          if(length(this_cp)>0){
               sigma_hat = sqrt(RSS/(2*n-2*m-2))
          }else{sigma_hat = sqrt(RSS/(2*n-4))}

          ## Record the velocity vector each segment and the
          if(length(this_cp)>0)
          {for (i in 1:m)
          {
               vx_vec[i] = sum(alpha[2:(i+1)])
          }
               vx_vec = as.vector(na.omit(vx_vec))

               for (i in 1:m)
               {
                    vy_vec[i] = sum(beta[2:(i+1)])
               }
               vy_vec = as.vector(na.omit(vy_vec))
               this_segvel = sqrt(vx_vec^2 + vy_vec^2)  # m elements
               # vel_vals = replicate(n,0)
               # if(m >2)
               # {vel_vals[1:this_cp[1]] = this_segvel[1]
               # for (i in 1:(m-2))
               # {vel_vals[this_cp[i]:this_cp[i+1]] = this_segvel[i+1]
               # }
               # vel_vals[this_cp[m-1]:length(t)] = this_segvel[m]
               #
               # ## time
               # this_segtime[1]=t[this_cp[1]]-t[1]
               # for (i in 2:(m-1))
               # {this_segtime[i]=t[this_cp[i]]-t[this_cp[i-1]]}
               # this_segtime[m]=t[n]-t[this_cp[m-1]]
               # }else{vel_vals[1:this_cp[1]] = this_segvel[1]
               # vel_vals[this_cp[1]:length(t)] = this_segvel[m]
               # this_segtime[1] = t[this_cp[1]]-t[1]
               # this_segtime[m] = t[n]-t[this_cp[1]]}

               this_segtime = diff(c(t[1],t[this_cp],
                                  t[length(t)]))


               ## eta=1/sigma_hat
               this_segeta = 1/(sigma_hat)

               ## theta
               for (i in 1:m)
               {
                    this_segtheta[i] = atan2(vy_vec[i],vx_vec[i])
               }
          }else{ ## For 0 changepoint
               vx_vec = alpha[2]
               vy_vec = beta[2]
               this_segvel = sqrt(vx_vec^2 + vy_vec^2)
               # vel_vals = c()
               # vel_vals[1:length(t)] = this_segvel
               this_segtime = t[n] - t[1]
               this_segeta = 1/(sigma_hat)
               this_segtheta = atan2(vy_vec,vx_vec)
          }

          if(old_version==TRUE){return(list(t = t,
                      x = x,
                      y = y,
                      x_piecewise = x_plinear,      # this is as long as the path
                      y_piecewise = y_plinear,      # this is as long as the path
                      this_cp = this_cp,
                      path_segvel = this_segvel,    # vector: one component per segment
                      path_segtime = this_segtime,  # vector: one component per segment
                      path_segtheta = this_segtheta,# vector: one component per segment
                      path_segeta = this_segeta,
                      inter_alpha = alpha_intercept,
                      inter_beta = beta_intercept,
                      u_x = vx_vec,
                      v_y = vy_vec,
                      path_RSS = RSS,
                      logical = TRUE
          ))}else{
               path_inferred = tibble(t = t,
                             j = rep("NA",length(t)),
                             x = x,
                             y = y,
                             a = x_plinear,
                             b = y_plinear)


               segments_inferred = tibble(
                         cp_times = c(t[this_cp],t[length(t)]),
                         durations = this_segtime,
                         states = infer_states_speed_cutoff(this_segvel,cutoff = 0.1),
                         speeds = this_segvel,
                         vx = vx_vec,
                         vy = vy_vec,
                         angles = this_segtheta)

               path_inferred$j = build_j_inferred(path_inferred,segments_inferred)

               return(list(
                    segments_inferred = segments_inferred,
                    path_inferred = path_inferred))

          }

     }

}


infer_states_speed_cutoff<-function(seg_speeds, cutoff = 0.1){
     return((seg_speeds > cutoff)*1)
}

build_j_inferred<-function(path_inf, seg_inf){
     j = 0*path_inf$t
     cpt = c(path_inf$t[1]-0.001,seg_inf$cp_times)
     for (i in 1:length(cpt)){
          seg_idx = which((path_inf$t > cpt[i])*(path_inf$t <= cpt[i+1]) == 1)
          j[seg_idx] = seg_inf$states[i]
     }
     return(j)
}

d_fn<-function(k,p)
{#k, p
     if (k<=p){d=0}else{d=1}
     return(d)
}



### Fill the missing data

fill_missing_data <- function(t, x, y, track_id)
{
     t = round(t,digits=4)
     this_dt = round(min(diff(t)),digits=4)
     diff_t = round(diff(t),digits=4)
     this_new_t = seq(t[1],t[length(t)],by=this_dt)
     this_new_x = 0*this_new_t
     this_new_y = 0*this_new_t

     sigma_est_x = sd(diff(x))/2
     sigma_est_y = sd(diff(y))/2

     find_next_good = function(x,idx) {
          i=idx
          while(x[i] == 0){
               i = i+1
          }
          return(i)
     }

     for (tt in 1:length(t)) {
          this_new_idx = which(abs((this_new_t - t[tt]))<0.0005)
          this_new_x[this_new_idx] = x[tt]
          this_new_y[this_new_idx] = y[tt]
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
          if ((tt > 1 && tt< length(this_new_t))|| tt < length(this_new_idx)) {
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
          if ((tt > 1 && tt< length(this_new_t)) || tt < length(this_new_idx)) {
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
     # plot(this_new_t,this_new_x,pch=20)
     # points(this_new_t[fill_idx_x],this_new_x[fill_idx_x],col="red",pch=20)
     # points(this_new_t[replace_idx_x],this_new_x[replace_idx_x],col="blue",pch=4)
     #
     # plot(this_new_t,this_new_y,pch=20)
     # points(this_new_t[fill_idx_y],this_new_y[fill_idx_y],col="red",pch=20)
     # points(this_new_t[replace_idx_y],this_new_y[replace_idx_y],col="blue",pch=4)

     return(list(t=this_new_t,
                 x=this_new_x,
                 y=this_new_y,
                 filled_x = fill_idx_x,
                 filled_y = fill_idx_y,
                 replaced_x = replace_idx_x,
                 replaced_y = replace_idx_y,
                 track_id = track_id))
}

