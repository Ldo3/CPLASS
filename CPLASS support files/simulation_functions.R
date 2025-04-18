get_segment = function(cur_state,cur_angle,theta){
  if (cur_state == 0){
    new_state = sample(c(0,1),1,prob = c(theta$p_SS,theta$p_SM))
  } else {
    new_state = sample(c(0,1),1,prob = c(theta$p_MS,theta$p_MM))
  }
  
  if (new_state == 0){
    new_angle = cur_angle
  } else {
    p_temp = runif(1)
    if (p_temp > (theta$p_reverse + theta$p_continue)){
      new_angle = runif(1,0,2*pi)
    } else if (p_temp > theta$p_reverse) {
      new_angle = cur_angle
    } else {
      new_angle = cur_angle + pi
      if(new_angle > 2*pi){
        new_angle = new_angle - 2*pi
      }
    }
  }
  new_sd = get_speed_duration(new_state,theta)

  return(tibble(
    state = new_state,
    angle = new_angle,
    speed = new_sd$speed,
    dur = new_sd$dur
  ))
}

get_speed_duration = function(cur_state, theta){
  if (cur_state == 0){
    new_speed = 0
    new_dur = rexp(1,rate = 1/theta$dur_stationary_avg)
  } else {
    new_speed = rgamma(1,theta$speed_alpha,rate = theta$speed_beta)
    if (theta$DEPENDENT_SPEED_DUR){
      new_dur = rexp(1, rate = new_speed/theta$dist_avg)
    } else {
      new_dur = rexp(1,rate = theta$speed_alpha/(theta_speed_beta*theta$dist_avg))
    }
  }
  return(tibble(
    speed = new_speed/1000,
    dur = new_dur
  ))
}

sim_segments = function(t,path_init,theta){
  t_final = t[length(t)]
  ##################
  # Determine the changepoint times and states
  # J(t) = 0 if stationary; J(t) = 1 if motile

  cp_times = t[1] + path_init$duration0
  seg_states = path_init$j0
  seg_durations = path_init$duration0
  seg_speeds = path_init$speed0
  seg_vx = path_init$speed0*cos(path_init$angle0)
  seg_vy = path_init$speed0*sin(path_init$angle0)
  seg_angles = path_init$angle0
  
  t_current = cp_times
  
  while(t_current < t_final){
    new_segment = get_segment(seg_states[length(seg_states)],
                              seg_angles[length(seg_angles)],
                              theta = theta)
    
    cp_times = c(cp_times,t_current + new_segment$dur)
    seg_states = c(seg_states,new_segment$state)
    seg_durations = c(seg_durations,new_segment$dur)
    seg_speeds = c(seg_speeds,new_segment$speed)
    seg_vx = c(seg_vx,new_segment$speed*cos(new_segment$angle))
    seg_vy = c(seg_vy,new_segment$speed*sin(new_segment$angle))
    seg_angles = c(seg_angles,new_segment$angle)

    t_current = t_current + new_segment$dur
  }
  
  num_states = length(cp_times)
  if (num_states == 1){
    seg_durations = t_final - t[1]
  }  else {
    seg_durations[length(seg_durations)] = t_final - cp_times[length(cp_times)-1]
  }
  
  seg_summary = tibble(
    cp_times = cp_times,
    durations = seg_durations,
    states = seg_states,
    speeds = seg_speeds,
    vx = seg_vx,
    vy = seg_vy,
    angles = seg_angles
  )

  return(seg_summary)
}

segments_to_path = function(t,path_init,seg_summary,theta){
  j = 0*t
  x = 0*t
  y = 0*t
  a = 0*t
  b = 0*t

  start_idx = 1
  seg_idx = which(t <= seg_summary$cp_times[start_idx])
  while(length(seg_idx) == 0){
    start_idx = start_idx + 1
    seg_idx = which(t <= seg_summary$cp_times[start_idx])
  }
  j[seg_idx] = rep(seg_summary$states[start_idx],length(seg_idx))

  ####### Need to add anchor diffusion #####
  a[seg_idx] = c(path_init$a0,path_init$a0 + 
    cumsum(seg_summary$vx[start_idx]*diff(t[seg_idx])))
  b[seg_idx] = c(path_init$b0,path_init$b0 + 
                   cumsum(seg_summary$vy[start_idx]*diff(t[seg_idx])))

  ###### Write observation points ######
  x[seg_idx] = a[seg_idx] + 
    theta$sigma_cargo*rnorm(length(seg_idx),0,1)
  y[seg_idx] = b[seg_idx] + 
    theta$sigma_cargo*rnorm(length(seg_idx),0,1)
  
  ###### Loop through segments #######
  a_last = a[seg_idx[length(seg_idx)]]
  b_last = b[seg_idx[length(seg_idx)]]
  for (i in (start_idx+1):length(seg_summary$cp_times)){
    seg_idx = which(t <= seg_summary$cp_times[i] & 
                      t > seg_summary$cp_times[i-1])
    if (length(seg_idx) > 0){
      j[seg_idx] = rep(seg_summary$states[i],length(seg_idx))
      
      t_seg = c(t[min(seg_idx)-1],t[seg_idx])
      a[seg_idx] = c(a_last + 
                       cumsum(seg_summary$vx[i]*diff(t_seg)))
      b[seg_idx] = c(b_last + 
                       cumsum(seg_summary$vy[i]*diff(t_seg)))
      x[seg_idx] = a[seg_idx] + 
        theta$sigma_cargo*rnorm(length(seg_idx),0,1)
      y[seg_idx] = b[seg_idx] + 
        theta$sigma_cargo*rnorm(length(seg_idx),0,1)
      
      a_last = a[max(seg_idx)]
      b_last = b[max(seg_idx)]
    }
  }

  path = tibble(
    t = t,
    j = j,
    x = x,
    y = y,
    a = a,
    b = b
  )
}

sim_path = function(Hz,n,theta,path_name = "Path",burnin_time = 20){
  Delta = 1/Hz
  
  t_final = Delta*n
  t = seq(0,t_final,by = Delta)
  t_burnin = seq(0,burnin_time,by = Delta)
  
  j0 = rbinom(1,1,theta$p_SM / (theta$p_SM + theta$p_MS))
  first_sd = get_speed_duration(j0,theta)
  burnin_init = tibble(
    a0 = 0,
    b0 = 0,
    j0 = j0,
    duration0 = first_sd$dur,
    speed0 = first_sd$speed,
    angle0 = runif(1)*2*pi
  )
  
  seg_burnin = sim_segments(t_burnin,burnin_init,theta)
  num_states = length(seg_burnin$states)
  
  path_init = tibble(
    a0 = 0,
    b0 = 0,
    j0 = seg_burnin$states[num_states],
    duration0 = seg_burnin$cp_times[num_states] - burnin_time,
    speed0 = seg_burnin$speeds[num_states],
    angle0 = seg_burnin$angles[num_states]
  )
  seg_summary = sim_segments(t,path_init,theta)
  this_path = segments_to_path(t,path_init,seg_summary,theta)
  path_data = list(
    segments = seg_summary,
    path = this_path,
    path_init = path_init,
    theta = theta,
    path_id = path_name
  )
  return(path_data)
}
