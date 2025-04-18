# Log-likelihood function
# library(extraDistr)
loglikelihood<- function(n,RSS)
{
     
     l=-n*log(2*pi)-n*log(RSS)+n*log(2*n)-n
     return(l)
}


CS<-function(t,x,y,r,p1,s_cap, with_speed)
{
     #p1: penalty for each parameter
     #s_cap: a minimum speed that has no penalty
     n=length(t)
     cps=which(r==1)+1
     # Compute log-likelihood
     pla= piecewise_linear_con(t,x,y,cps)
     if(pla$logical==TRUE){
          llh= -n*log(2*pi)-n*log(pla$path_RSS)+n*log(2*n)-n
               # -n*log(pla$path_RSS)
          seg_vel = pla$path_segvel
      
          # pv = 0
          # for(i in 1: length(seg_vel))
          # {
          #      pv = pv + max(1,exp(nu*(seg_vel[i]-s_cap)))
          # }
          
          pv = 0
          for(i in 1: length(seg_vel))
          {
               pv = pv + max(0,(seg_vel[i]-s_cap))
          }

          p = p1 * length(cps) + with_speed* pv 
          l = llh - p
          # k=length(cps)
          # pen(mo3$d)-1.25*(k-length(mo3$m3))
          # p=pen(mo3$d)+(k-length(mo3$m3))
          return(list(s=l,
                      llh=llh,
                      p=p,
                      pv = pv,
                      logical=TRUE))
     }else{return(list(logical=FALSE))}
     
}

cp_to_r = function(cp,n){
     output = 0*(1:(n-2))
     output[cp-1] = 1
     return(output)
}

r_to_cp = function(r){
     1 + which(r == 1)
}


