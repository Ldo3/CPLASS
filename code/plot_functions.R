col_actual = "orange"
col_inferred = "steelblue"
col_stationary = "red3"
col_motile = "green4"
col_stationary_inferred = "pink"
col_motile_inferred = "green"

get_one_plot_size = function(path_info){
  path = path_info$path_inferred
  dist_max = max(c(max(path$x)-min(path$x)),
                            max((path$y)-min(path$y)))
  t_max = max(path$t) - min(path$t)

  return(tibble(
    xy_width = 1.2*max(dist_max),
    t_begin = 0,
    t_end = t_max
  ))
}

get_plot_size = function(path_list){
  dist_max = c()
  t_max = c()
  for (i in 1:length(path_list)){
    path = path_list[[i]]$path
    dist_max = c(dist_max,max(c(max(path$x)-min(path$x)),
                              max((path$y)-min(path$y))))
    t_max = c(t_max,max(path$t) - min(path$t))
  }

  return(tibble(
    xy_width = 1.2*max(dist_max),
    t_begin = 0,
    t_end = max(t_max)
  ))
}

get_plot_size_inferred = function(path_list){
  dist_max = c()
  t_max = c()
  for (i in 1:length(path_list)){
    path = path_list[[i]]$path_inferred
    dist_max = c(dist_max,max(c(max(path$x)-min(path$x)),
                              max((path$y)-min(path$y))))
    t_max = c(t_max,max(path$t) - min(path$t))
  }

  return(tibble(
    xy_width = 1.2*max(dist_max),
    t_begin = 0,
    t_end = max(t_max)
  ))
}


plot_path_inferred <- function(path_info,xy_width = NA,t_lim = NA, motor, max_speed, state_shaded = TRUE, title_ind = TRUE, show_time_changes = FALSE){
     if (is.na(xy_width) | is.na(t_lim[1])){
          frame_info = get_plot_size(list(path_info))
     }
     if (is.na(xy_width)){
          xy_width = frame_info$xy_width
     }
     if (is.na(t_lim[1])){
          t_lim = c(frame_info$t_begin,frame_info$t_end)
     }

     path = path_info$path_inferred
     segments = path_info$segments_inferred
     x_min = mean(path$x)-xy_width/2
     x_max = mean(path$x)+xy_width/2
     y_min = mean(path$y)-xy_width/2
     y_max = mean(path$y)+xy_width/2

     #################
     # This plot has x vs y and a vs b
     if(title_ind){p_xy = ggplot() +
          geom_path(data = path,
                    aes(x = x, y = y))+
          geom_path(data = path,
                    aes(x = a, y = b), col=col_inferred) +
          xlim(x_min,x_max)+
          ylim(y_min,y_max)+
          ggtitle(paste0(motor, " Path ",i)) +
          #for real data
          #ggtitle(paste0(motor, " Path ",i)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))}else{
               p_xy = ggplot() +
                    geom_path(data = path,
                              aes(x = x, y = y))+
                    geom_path(data = path,
                              aes(x = a, y = b), col=col_inferred) +
                    xlim(x_min,x_max)+
                    ylim(y_min,y_max)+
                    ggtitle(paste0(motor, " Path ")) +
                    #for real data
                    #ggtitle(paste0(motor, " Path ",i)) +
                    theme_minimal() +
                    theme(plot.title = element_text(hjust = 0.5))
          }
     

     # This creates the (x,a) vs t and (y,b) vs t plots
     p_xt = ggplot()
     p_yt = ggplot()
     if (state_shaded){
          shades = c(col_stationary_inferred,col_motile_inferred)
          #make sure the cp_time is not longer than the path length
          if(segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]){
               segments$cp_times[length(segments$cp_times)] = path$t[length(path$t)]
          }
          seg_ends = c(path$t[1],segments$cp_times)
          for (i in 1:length(segments$cp_times)) {
               p_xt = p_xt + annotate("rect",
                                      xmin = seg_ends[i], xmax = seg_ends[i+1],
                                      ymin = x_min, ymax = x_max,
                                      alpha = 0.2, fill = shades[segments$states[i]+1]
               )

               p_yt = p_yt + annotate("rect",
                                      xmin = seg_ends[i], xmax = seg_ends[i+1],
                                      ymin = y_min, ymax = y_max,
                                      alpha = 0.2, fill = shades[segments$states[i]+1]
               )
          }
     }

     if(show_time_changes == FALSE){
     p_xt = p_xt +
          #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
          #                 col="pink",linewidth = 0.25) +
          geom_path(data = path,
                    aes(x = t, y = x))+
          geom_path(data = path,
                    aes(x = t, y = a), col=col_inferred) +
          xlim(t_lim[1],t_lim[2]) +
          ylim(x_min,x_max) +
          ggtitle("x- and a- Time Series") +
          theme_minimal() +
          theme(axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                plot.title = element_text(hjust = 0.5))

     p_yt = p_yt +
          #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
          #                 col="pink",linewidth = 0.25) +
          geom_path(data = path,
                    aes(x = t, y = y))+
          geom_path(data = path,
                    aes(x = t, y = b), col=col_inferred) +
          xlim(t_lim[1],t_lim[2]) +
          ylim(y_min,y_max) +
          ggtitle("y- and b- Time Series") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
     }else{
          p_xt = p_xt +
                    geom_vline(xintercept = segments$cp_times[-length(segments$cp_times)], linetype="twodash",
                               col="black",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = x))+
               geom_path(data = path,
                         aes(x = t, y = a), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(x_min,x_max) +
               ggtitle("x- and a- Time Series") +
               theme_minimal() +
               theme(axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(hjust = 0.5))
          
          p_yt = p_yt +
                    geom_vline(xintercept = segments$cp_times[-length(segments$cp_times)], linetype="twodash",
                               col="black",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = y))+
               geom_path(data = path,
                         aes(x = t, y = b), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(y_min,y_max) +
               ggtitle("y- and b- Time Series") +
               theme_minimal() +
               theme(plot.title = element_text(hjust = 0.5))
     }

     p_segments = ggplot() +
          geom_hline(yintercept = 0, linewidth = 0.5, col="gray") +
          geom_vline(xintercept = 0, linewidth = 0.5, col="gray") +
          geom_point(data = segments, aes(x = durations, y = speeds),
                     alpha = 0.4,size = 0.6,col = col_inferred) +
          xlim(0,t_lim[2]) + ylim(0,max_speed) + xlab("Segment duration (s)") +
          ylab(expression(paste("Speed (", mu, "m/s)")))+
          ggtitle(paste0("Segments: ",length(segments$cp_times))) +
          theme_classic()+
          theme(plot.title = element_text(hjust = 1, size = 10,
                                          margin = margin(0,0,-10,0)))

     layout <- c(
          area(t = 1, l = 1, b = 5, r = 3),
          area(t = 1, l = 4, b = 4, r = 7),
          area(t = 5, l = 4, b = 7, r = 7),
          area(t = 6, l = 1, b = 7, r = 3)
     )

     path_dashboard = (p_xy + p_xt + p_yt + p_segments +
                            plot_layout(design = layout))

     return(path_dashboard)
}

plot_path_inferred_second_version <- function(path_info,xy_width = NA,t_lim = NA, motor, max_speed, state_shaded = TRUE, title_ind = TRUE, show_time_changes = FALSE){
     if (is.na(xy_width) | is.na(t_lim[1])){
          frame_info = get_plot_size(list(path_info))
     }
     if (is.na(xy_width)){
          xy_width = frame_info$xy_width
     }
     if (is.na(t_lim[1])){
          t_lim = c(frame_info$t_begin,frame_info$t_end)
     }
     
     path = path_info$path_inferred
     segments = path_info$segments_inferred
     # x_min = mean(path$x)-xy_width 
     x_min = min(path$x) -0.01
     # x_max = mean(path$x)+xy_width
     x_max = max(path$x) +0.01
     # y_min = mean(path$y)-xy_width
     # y_max = mean(path$y)+xy_width
     y_min = min(path$y) -0.01
     y_max = max(path$y) +0.01
     
     #################
     # This plot has x vs y and a vs b
     if(title_ind){p_xy = ggplot() +
          geom_path(data = path,
                    aes(x = x, y = y))+
          geom_path(data = path,
                    aes(x = a, y = b), col=col_inferred) +
          xlim(x_min,x_max)+
          ylim(y_min,y_max)+
          ggtitle(paste0(motor, " Path ",i)) +
          #for real data
          #ggtitle(paste0(motor, " Path ",i)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))}else{
               p_xy = ggplot() +
                    geom_path(data = path,
                              aes(x = x, y = y))+
                    geom_path(data = path,
                              aes(x = a, y = b), col=col_inferred) +
                    xlim(x_min,x_max)+
                    ylim(y_min,y_max)+
                    ggtitle(paste0(motor, " Path ")) +
                    #for real data
                    #ggtitle(paste0(motor, " Path ",i)) +
                    theme_minimal() +
                    theme(plot.title = element_text(hjust = 0.5))
          }
     
     
     # This creates the (x,a) vs t and (y,b) vs t plots
     p_xt = ggplot()
     p_yt = ggplot()
     if (state_shaded){
          shades = c(col_stationary_inferred,col_motile_inferred)
          #make sure the cp_time is not longer than the path length
          if(segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]){
               segments$cp_times[length(segments$cp_times)] = path$t[length(path$t)]
          }
          seg_ends = c(path$t[1],segments$cp_times)
          for (i in 1:length(segments$cp_times)) {
               p_xt = p_xt + annotate("rect",
                                      xmin = seg_ends[i], xmax = seg_ends[i+1],
                                      ymin = x_min, ymax = x_max,
                                      alpha = 0.2, fill = shades[segments$states[i]+1]
               )
               
               p_yt = p_yt + annotate("rect",
                                      xmin = seg_ends[i], xmax = seg_ends[i+1],
                                      ymin = y_min, ymax = y_max,
                                      alpha = 0.2, fill = shades[segments$states[i]+1]
               )
          }
     }
     
     if(show_time_changes == FALSE){
          p_xt = p_xt +
               #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
               #                 col="pink",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = x))+
               geom_path(data = path,
                         aes(x = t, y = a), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(x_min,x_max) +
               ggtitle("x- and a- Time Series") +
               theme_minimal() +
               theme(axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(hjust = 0.5))
          
          p_yt = p_yt +
               #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
               #                 col="pink",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = y))+
               geom_path(data = path,
                         aes(x = t, y = b), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(y_min,y_max) +
               ggtitle("y- and b- Time Series") +
               theme_minimal() +
               theme(plot.title = element_text(hjust = 0.5))
     }else{
          p_xt = p_xt +
               geom_vline(xintercept = segments$cp_times[-length(segments$cp_times)], linetype="twodash",
                          col="black",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = x))+
               geom_path(data = path,
                         aes(x = t, y = a), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(x_min,x_max) +
               ggtitle("x- and a- Time Series") +
               theme_minimal() +
               theme(axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(hjust = 0.5))
          
          p_yt = p_yt +
               geom_vline(xintercept = segments$cp_times[-length(segments$cp_times)], linetype="twodash",
                          col="black",linewidth = 0.25) +
               geom_path(data = path,
                         aes(x = t, y = y))+
               geom_path(data = path,
                         aes(x = t, y = b), col=col_inferred) +
               xlim(t_lim[1],t_lim[2]) +
               ylim(y_min,y_max) +
               ggtitle("y- and b- Time Series") +
               theme_minimal() +
               theme(plot.title = element_text(hjust = 0.5))
     }
     
     p_segments = ggplot() +
          geom_hline(yintercept = 0, linewidth = 0.5, col="gray") +
          geom_vline(xintercept = 0, linewidth = 0.5, col="gray") +
          geom_point(data = segments, aes(x = durations, y = speeds),
                     alpha = 0.4,size = 0.6,col = col_inferred) +
          xlim(0,t_lim[2]) + ylim(0,max_speed) + xlab("Segment duration (s)") +
          ylab(expression(paste("Speed (", mu, "m/s)")))+
          ggtitle(paste0("Segments: ",length(segments$cp_times))) +
          theme_classic()+
          theme(plot.title = element_text(hjust = 1, size = 10,
                                          margin = margin(0,0,-10,0)))
     
     layout <- c(
          area(t = 1, l = 1, b = 1, r = 2),  # p_xy in row 1, col 1–2
          area(t = 1, l = 3, b = 1, r = 4),  # p_segments in row 1, col 3–4
          area(t = 2, l = 1, b = 2, r = 4),  # p_xt in row 2
          area(t = 3, l = 1, b = 3, r = 4)   # p_yt in row 3
     )
     
     path_dashboard <- (
          p_xy + p_segments + p_xt + p_yt +
               plot_layout(design = layout)
     )
     
     return(path_dashboard)
}

plot_path_inferred_xy <- function(path_info,xy_width = NA,t_lim = NA, motor, max_speed, state_shaded = TRUE){
     if (is.na(xy_width) | is.na(t_lim[1])){
          frame_info = get_plot_size(list(path_info))
     }
     if (is.na(xy_width)){
          xy_width = frame_info$xy_width
     }
     if (is.na(t_lim[1])){
          t_lim = c(frame_info$t_begin,frame_info$t_end)
     }

     path = path_info$path_inferred
     segments = path_info$segments_inferred
     x_min = mean(path$x)-xy_width/2
     x_max = mean(path$x)+xy_width/2
     y_min = mean(path$y)-xy_width/2
     y_max = mean(path$y)+xy_width/2

     #################
     # This plot has x vs y and a vs b
     p_xy = ggplot() +
          geom_path(data = path,
                    aes(x = x, y = y))+
          geom_path(data = path,
                    aes(x = a, y = b), col=col_inferred) +
          xlim(x_min,x_max)+
          ylim(y_min,y_max)+
          ggtitle(paste0(motor, " Path ")) +

          #for real data
          #ggtitle(paste0(motor, " Path ",i)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))

     # This creates the (x,a) vs t and (y,b) vs t plots

     path_dashboard = p_xy

     return(path_dashboard)
}

plot_path_actual_and_inferred = function(path_info,xy_width = NA,
                                         t_lim = NA, state_shaded = TRUE){
  if (is.na(xy_width) | is.na(t_lim[1])){
    frame_info = get_plot_size(path_info)
  }
  if (is.na(xy_width)){
    xy_width = frame_info$xy_width
  }
  if (is.na(t_lim[1])){
    t_lim = c(frame_info$t_begin,frame_info$t_end)
  }

  #actual
  path = path_info$path
  segments = path_info$segments

  #inferred
  path_inf = path_info$path_inferred
  segments_inf = path_info$segments_inferred

  duration_difference = duration_of_difference(path_info)

  x_min = mean(path$x)-xy_width/2
  x_max = mean(path$x)+xy_width/2
  y_min = mean(path$y)-xy_width/2
  y_max = mean(path$y)+xy_width/2

  #################
  # This plot has x vs y and a vs b
  p_xy = ggplot() +
    geom_path(data = path,
              aes(x = x, y = y))+
    geom_path(data = path,
              aes(x = a, y = b), col=col_actual) +
    geom_path(data = path_inf,
              aes(x = a, y = b), col=col_inferred) +
    xlim(x_min,x_max)+
    ylim(y_min,y_max)+
    ggtitle(paste0("Simulated path ",i)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # This creates the (x,a) vs t and (y,b) vs t plots
  p_xt = ggplot()
  p_yt = ggplot()
  if (state_shaded){
    #actual
    shades = c(col_stationary,col_motile)
    #make sure the cp_time is not longer than the path length
    if(segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]){
      segments$cp_times[length(segments$cp_times)] = path$t[length(path$t)]
    }
    seg_ends = c(path$t[1],segments$cp_times)
    for (i in 1:length(segments$cp_times)) {
      p_xt = p_xt + annotate("rect",
                             xmin = seg_ends[i], xmax = seg_ends[i+1],
                             ymin = x_min, ymax = x_max,
                             alpha = 0.2, fill = shades[segments$states[i]+1]
      )

      p_yt = p_yt + annotate("rect",
                             xmin = seg_ends[i], xmax = seg_ends[i+1],
                             ymin = y_min, ymax = y_max,
                             alpha = 0.2, fill = shades[segments$states[i]+1]
      )
    }

    #inferred
    shades = c(col_stationary_inferred,col_motile_inferred)
    #make sure the cp_time is not longer than the path length
    if(segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]){
      segments$cp_times[length(segments$cp_times)] = path$t[length(path$t)]
    }
    seg_ends = c(path_inf$t[1],segments_inf$cp_times)
    for (i in 1:length(segments_inf$cp_times)) {
      p_xt = p_xt + annotate("rect",
                             xmin = seg_ends[i], xmax = seg_ends[i+1],
                             ymin = x_min, ymax = x_max,
                             alpha = 0.2, fill = shades[segments_inf$states[i]+1]
      )

      p_yt = p_yt + annotate("rect",
                             xmin = seg_ends[i], xmax = seg_ends[i+1],
                             ymin = y_min, ymax = y_max,
                             alpha = 0.2, fill = shades[segments_inf$states[i]+1]
      )
    }
  }

  p_xt = p_xt +
    #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
    #                 col="pink",linewidth = 0.25) +
    geom_path(data = path,
              aes(x = t, y = x))+
    geom_path(data = path,
              aes(x = t, y = a), col=col_actual) +
    geom_path(data = path_inf,
              aes(x = t, y = a), col=col_inferred) +
    xlim(t_lim[1],t_lim[2]) +
    ylim(x_min,x_max) +
    ggtitle("x- and a- Time Series") +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))

  p_yt = p_yt +
    #      geom_vline(xintercept = segments$cp_times, linetype="twodash",
    #                 col="pink",linewidth = 0.25) +
    geom_path(data = path,
              aes(x = t, y = y))+
    geom_path(data = path,
              aes(x = t, y = b), col=col_actual) +
    geom_path(data = path_inf,
              aes(x = t, y = b), col=col_inferred) +
    xlim(t_lim[1],t_lim[2]) +
    ylim(y_min,y_max) +
    ggtitle("y- and b- Time Series") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  p_segments = ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.5, col="gray") +
    geom_vline(xintercept = 0, linewidth = 0.5, col="gray") +
    geom_point(data = segments, aes(x = durations, y = speeds),
               alpha = 0.4,size = 0.6,col = col_actual) +
    geom_point(data = segments_inf, aes(x = durations, y = speeds),
               alpha = 0.4,size = 0.6,col = col_inferred) +
    xlim(t_lim[1],t_lim[2]) + ylim(0,1) + xlab("Segment duration (s)") +
       ylab(expression(paste("Speed (", mu, "m/s)")))+
    ggtitle(paste0("Segments: ",length(segments$cp_times),
                "\n Inf. Segments: ",length(segments_inf$cp_times),
                "\n inference Gap: ",round(duration_difference,digits=2))) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 1, size = 10,
                                    margin = margin(0,0,-10,0)))

  layout <- c(
    area(t = 1, l = 1, b = 5, r = 3),
    area(t = 1, l = 4, b = 4, r = 7),
    area(t = 5, l = 4, b = 7, r = 7),
    area(t = 6, l = 1, b = 7, r = 3)
  )

  path_dashboard = (p_xy + p_xt + p_yt + p_segments +
                      plot_layout(design = layout))

  return(path_dashboard)
}

plot_csa = function(csa, max_error = 0.1, legend = TRUE){
  Hz = csa %>% distinct(Hz)
  pl_title = paste("Cumulative Speed Allocation:", csa$Hz[which(!is.na(csa$Hz))],"Hz")
  p_csa = ggplot()
  if (legend == TRUE){
    csa_three = csa %>% filter(label == "True" & s > 0)
    csa_three = bind_rows(csa_three, csa %>% filter(s > 0 & subsample == 1))
    p_csa = p_csa +
      geom_path(data = csa_three,
                aes(x = s, y = csa, group = label, col = label))
  }
  p_csa = p_csa +
    geom_path(data = csa %>% filter(label == "Sample" & s > 0),
              aes(x = s, y = csa, group = subsample), col=col_actual,
              alpha = 0.25)+
    geom_path(data = csa %>% filter(label == "Inferred" & s > 0),
              aes(x = s, y = csa, group = subsample), col=col_inferred,
              alpha = 0.25)+
    geom_path(data = csa %>% filter(label == "True" & s > 0),
              aes(x = s, y = csa), col="black")+
    ggtitle(pl_title)+
    coord_cartesian(ylim = c(0,1))+
    theme_classic()+
    theme(legend.position = c(0.9,0.2))

  return(p_csa)
}


# Loop through the list and add rectangles to the plot
#for (rect in rectangles) {
#  plot <- plot +
#    annotate(
#      "rect",
#      xmin = rect$xmin, xmax = rect$xmax,
#      ymin = rect$ymin, ymax = rect$ymax,
#      alpha = 0.2, fill = rect$fill
#    )
#}


duration_of_difference = function(path_info){
     cp_times = path_info$segments$cp_times
     cp_times[length(cp_times)] = max(path_info$path$t)
     segments1 <- data.frame(start = c(path_info$path$t[1],
                                       cp_times[-length(cp_times)]), 
                             end = cp_times,
                             state = path_info$segments$states)
     cp_times = path_info$segments_inferred$cp_times
     segments2 <- data.frame(start = c(path_info$path$t[1],
                                       cp_times[-length(cp_times)]), 
                             end = cp_times,
                             state = path_info$segments_inferred$states)
     
     
     total_overlap_length <- 0
     
     i = 1
     j = 1
     
     while((i <= length(segments1[,1])) && (j <= length(segments2[,1]))){ 
          
          count1 = 1
          count2 = 1
          chain1 = segments1[i,]
          chain2 = segments2[j,]
          
          if(chain1$end == chain2$end){
               if(chain1$state != chain2$state){
                    # Calculate the overlap length
                    overlap_length <- min(chain1$end, chain2$end) - max(chain1$start,chain2$start)
                    
                    # Add the length to the total overlap
                    total_overlap_length <- total_overlap_length + overlap_length
               }
               i = i + 1
               j = j + 1
          }
          
          if(chain1$end < chain2$end){
               if(chain1$state != chain2$state){
                    # Calculate the overlap length
                    overlap_length <- min(chain1$end, chain2$end) - max(chain1$start,chain2$start)
                    
                    # Add the length to the total overlap
                    total_overlap_length <- total_overlap_length + overlap_length
               }
               #check to see if next segment in chain1 is still less than chain2
               next_chain1 = segments1[i+count1,]
               while(next_chain1$end <= chain2$end){
                    if(next_chain1$state != chain2$state){
                         # Calculate the overlap length
                         overlap_length <- next_chain1$end - next_chain1$start
                         
                         # Add the length to the total overlap
                         total_overlap_length <- total_overlap_length + overlap_length
                    }
                    count1 = count1 + 1
                    if((i+count1) > length(segments1[,1])){
                         break
                    }
                    next_chain1 = segments1[i+count1,]
               }
               if((i+count1) > length(segments1[,1])){
                    break
               }  else if(next_chain1$start <= chain2$end){
                    if(next_chain1$state != chain2$state){
                         # Calculate the overlap length
                         overlap_length <- chain2$end - next_chain1$start
                         
                         # Add the length to the total overlap
                         total_overlap_length <- total_overlap_length + overlap_length
                    }
               }
               i = i + count1
               j = j + 1
               
          }
          
          if(chain1$end > chain2$end){
               if(chain1$state != chain2$state){
                    # Calculate the overlap length
                    overlap_length <- min(chain1$end, chain2$end) - max(chain1$start,chain2$start)
                    
                    # Add the length to the total overlap
                    total_overlap_length <- total_overlap_length + overlap_length
               }
               #check to see if next segment in chain2 is still less than chain1
               next_chain2 = segments2[j+count2,]
               while(next_chain2$end <= chain1$end){
                    if(next_chain2$state != chain1$state){
                         # Calculate the overlap length
                         overlap_length <- next_chain2$end - next_chain2$start
                         
                         # Add the length to the total overlap
                         total_overlap_length <- total_overlap_length + overlap_length
                    }
                    count2 = count2 + 1
                    if((j+count2) > length(segments2[,1])){
                         break
                    }
                    next_chain2 = segments2[j+count2,]
               }    
               if((j+count2) > length(segments2[,1])){
                    break
               }  else if(next_chain2$start <= chain1$end){
                    if(next_chain2$state != chain1$state){
                         # Calculate the overlap length
                         overlap_length <- chain1$end - next_chain2$start
                         
                         # Add the length to the total overlap
                         total_overlap_length <- total_overlap_length + overlap_length
                    }
               }
               j = j + count2
               i = i + 1
               
          }
          
     }
     
     return(total_overlap_length)
}

