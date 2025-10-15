# Log-likelihood function
# library(extraDistr)
loglikelihood<- function(n,RSS, sd = NA)
{
     if(!is.na(sd)){
       l = -n*log*(2*pi) - n*log(sd^2) - 1/(2*sd^2) - 1/(2*sd^2)*RSS
     }else{
       l=-n*log(2*pi)-n*log(RSS)+n*log(2*n)-n
     }

     return(l)
}

AICc = function(ncp, n)
{
  k = 2*(ncp+2)+1
  if(2*n-k-1<0){aicc = 0}else{
  aicc = 2*k + 2*k*(k+1)/(2*n-k-1)}
  return(aicc)
}

sSIC = function(ncp, n, gamma)
{

  k = 2*(ncp+2)+1
  ssic = log(n)^{gamma}*k
  return(ssic)
}

CS<-function(t,x,y,r,p1,s_cap, with_speed, sd = NA, pen="ssic", gamma)
{
     #p1: penalty for each parameter
     #s_cap: a minimum speed that has no penalty
     n=length(t)
     cps=which(r==1)+1
     # Compute log-likelihood
     pla= piecewise_linear_con(t,x,y,cps)
     if(pla$logical==TRUE){
          # llh= -n*log(2*pi)-n*log(pla$path_RSS)+n*log(2*n)-n
          llh = loglikelihood(n, pla$path_RSS, sd= sd)
               # -n*log(pla$path_RSS)
          seg_vel = pla$path_segvel

          pv = 0
          for(i in 1: length(seg_vel))
          {
               pv = pv + max(0,(seg_vel[i]-s_cap))
          }

          ncp = length(cps)

          if(pen == "hybrid"){
            p = max(AICc(ncp,n), sSIC(ncp,n, gamma)) + with_speed* pv
            }else if(pen == "aicc"){
                  p = AICc(ncp,n) + with_speed* pv
                  }else if(pen == "ssic"){
                       p = sSIC(ncp,n, gamma) + with_speed* pv}
          l = llh - p
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


