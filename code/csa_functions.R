compute_csa = function(segments_summary,speed_mesh){
  wecdf = 0*speed_mesh
  total_time = sum(segments_summary$durations)
  for (i in 1:length(speed_mesh)){
    wecdf[i] = sum(segments_summary$durations*
                     (segments_summary$speeds < speed_mesh[i]))/total_time
  }
  return(tibble(s = speed_mesh,
                csa = wecdf))
}

csa_theoretical = function(theta, s_mesh){
  alpha = theta$speed_alpha
  beta = theta$speed_beta
  c = 1/theta$dist_avg
  p_MM = theta$p_MM
  dist_avg = theta$dist_avg
  num_M_avg = (1/(1-p_MM))

  dur_S_avg = theta$dur_stationary_avg
  if (theta$DEPENDENT_SPEED_DUR){
    dur_M_avg = beta/(c*(alpha-1))
    exp_delta_I = dur_S_avg + num_M_avg*dur_M_avg * pgamma(s_mesh*1000,shape = alpha-1, rate = beta)
  } else{
    dur_M_avg = dist_avg/(alpha/beta)
    exp_delta_I = dur_S_avg + num_M_avg*dur_M_avg * pgamma(s_mesh*1000,shape = alpha, rate = beta)
  }

  exp_delta_T = dur_S_avg + num_M_avg*dur_M_avg

  return(exp_delta_I/exp_delta_T)
}

summarize_segments  = function(path_list,group_label){
  segments_summary = tibble(
    durations = c(),
    states = c(),
    speeds = c(),
    label = c(),
    path_id = c()
  )

  for (i in 1:length(path_list)){
    this_path = path_list[[i]]

    these_segments = this_path$segments %>% select(durations, states, speeds)
    num_segments = length(these_segments$durations)

    these_segments = these_segments %>%
      mutate(label = rep(group_label,length(this_path$segments$states)),
             path_id = rep(i,length(this_path$segments$states)))
    segments_summary = bind_rows(segments_summary, these_segments)
  }
  return(segments_summary)
}

summarize_segments_inferred  = function(path_list,group_label){
  segments_summary = tibble(
    durations = c(),
    states = c(),
    speeds = c(),
    label = c(),
    path_id = c()
  )

  for (i in 1:length(path_list)){
    this_path = path_list[[i]]

    these_segments = this_path$segments_inferred %>% select(durations, states, speeds)
    num_segments = length(these_segments$durations)

    these_segments = these_segments %>%
      mutate(label = rep(group_label,length(this_path$segments_inferred$states)),
             path_id = rep(i,length(this_path$segments_inferred$states)))
    segments_summary = bind_rows(segments_summary, these_segments)
  }
  return(segments_summary)
}
summarize_segments_test  = function(path_list,group_label){
  segments_summary = tibble(
    durations = c(),
    states = c(),
    speeds = c(),
    label = c(),
    path_id = c()
  )

  for (i in 1:length(path_list)){
    this_path = path_list[[i]]

    these_segments = this_path$segments_inferred %>% select(durations, states_test, speeds)
    num_segments = length(these_segments$durations)

    these_segments = these_segments %>%
      mutate(label = rep(group_label,num_segments),
             path_id = rep(i,num_segments))
    these_segments = these_segments %>% rename(states = states_test)
    segments_summary = bind_rows(segments_summary, these_segments)
  }
  return(segments_summary)
}

summarize_segments_cutoff  = function(path_list,group_label){
  segments_summary = tibble(
    durations = c(),
    states = c(),
    speeds = c(),
    label = c(),
    path_id = c()
  )

  for (i in 1:length(path_list)){
    this_path = path_list[[i]]

    these_segments = this_path$segments_inferred %>% select(durations, states, speeds)
    num_segments = length(these_segments$durations)

    these_segments = these_segments %>%
      mutate(label = rep(group_label,length(this_path$segments_inferred$states)),
             path_id = rep(i,length(this_path$segments_inferred$states)))
    segments_summary = bind_rows(segments_summary, these_segments)
  }
  return(segments_summary)
}
