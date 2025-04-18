# source("piecewise_linear2.R")
source("CPLASS.R")
library("extraDistr") # this is for the Gumbel distribution

data = readRDS("CPLASS_21PN.rds")
pl = list()
speed_pen = c("without", "with")
with_speed = 0 # Choose 0 for without speed penalty, Choose 1 for with speed penalty

pen_name = c("2+lnn", "Gumbel", "4")
pen_coefficient = 1 # Choose 0 for 2+log(n), Choose 1 for Gumbel (g + log(n/10)), Choose 2 for old penalty (4)

print(paste("WE ARE USING THE", pen_name[pen_coefficient+1], "COEFFICIENT", speed_pen[with_speed+1], "SPEED PENALTY."))


i=1
while (i <= length(data))
{ #i=2
  print(paste("we are in the case ", i, sep = ""))
  # this_t = data[[i]]$this_t
  # this_x = data[[i]]$x_sim[[1]]
  # this_y = data[[i]]$y_sim[[1]]
  n = length(data[[i]]$this_t)
  
  ### Choose penalty coefficient p1
  pen_value = c(2+log(n), qgumbel(0.95) + log(n/10), 4 )
  p1 = pen_value[pen_coefficient+1]
  
  ###
  this_t = data[[i]]$this_t[-c(1, n)] #delete the starting and ending points of the paths
  this_x = data[[i]]$this_x[-c(1, n)]
  this_y = data[[i]]$this_y[-c(1, n)]
  time_rate = 0.05
  MCMC = try(MHpen(this_t,this_x,this_y,time_rate = time_rate,iter_max = 5000, burn_in = 500, p1=p1, v0=1, with_speed = with_speed), silent = FALSE)
  if(inherits(MCMC, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    #next
    i=i
  }else{
    Final_MCMC= FinalMH(MCMC)
    index.max = which.max(Final_MCMC$slc)
    pl[[i]]= piecewise_linear_con(this_t, this_x, this_y, Final_MCMC$this_cp[[index.max]])
    i = i+1
  }
  # if((i%%100)==0){outfile = paste("yap_list21PN_",i,sep="")
  # saveRDS(pl,file= outfile)}
  
}


outfile = paste("CPLASS_21PN_pen_", pen_name[pen_coefficient+1],"_", speed_pen[with_speed+1],"_speed_pen.rds",sep="")
saveRDS(pl,file= outfile)
# 

######