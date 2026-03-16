peak_ranges <- function(df, min_height = 0.25, min_time_distance = 0.1, expand_frac = 0.05) {
  library(pracma)
  library(dplyr)
  library(tibble)
  
  df <- df |> arrange(time)
  dt <- median(diff(df$time))
  min_distance_rows <- round(min_time_distance / dt)
  
  # Find peaks
  p <- findpeaks(df$intensity, minpeakheight = min_height, minpeakdistance = min_distance_rows)
  if (is.null(p)) return(tibble())
  
  apex_indices <- p[,2]
  
  # Sort peaks by apex
  apex_indices <- sort(apex_indices)
  
  # Expand peaks but prevent overlap
  peak_ranges_list <- vector("list", length(apex_indices))
  
  for(i in seq_along(apex_indices)) {
    apex_idx <- apex_indices[i]
    peak_val <- df$intensity[apex_idx]
    threshold <- peak_val * expand_frac
    
    # expand backward
    start_idx <- apex_idx
    while(start_idx > 1 && df$intensity[start_idx] > threshold) start_idx <- start_idx - 1
    # don't go past previous apex
    if(i > 1) start_idx <- max(start_idx, apex_indices[i-1] + 1)
    
    # expand forward
    end_idx <- apex_idx
    while(end_idx < nrow(df) && df$intensity[end_idx] > threshold) end_idx <- end_idx + 1
    # don't go past next apex
    if(i < length(apex_indices)) end_idx <- min(end_idx, apex_indices[i+1] - 1)
    
    peak_ranges_list[[i]] <- c(start_idx, end_idx)
  }
  
  peak_ranges_mat <- do.call(rbind, peak_ranges_list)
  
  t <- tibble(
    start_time  = df$time[peak_ranges_mat[,1]],
    apex_time   = df$time[apex_indices],
    end_time    = df$time[peak_ranges_mat[,2]],
    peak_height = p[order(apex_indices),1]
  )
  
  # Keep all peaks separate; no merging unless apexes are extremely close
  if(nrow(t) == 0) return(t)
  
  merged <- list(t[1,])
  
  for(i in 2:nrow(t)) {
    prev <- merged[[length(merged)]]
    apex_distance <- t$apex_time[i] - prev$apex_time
    
    if(apex_distance < 0.02) {  # tiny bumps merge
      merged[[length(merged)]] <- tibble(
        start_time  = min(prev$start_time, t$start_time[i]),
        apex_time   = if(t$peak_height[i] > prev$peak_height) t$apex_time[i] else prev$apex_time,
        end_time    = max(prev$end_time, t$end_time[i]),
        peak_height = max(prev$peak_height, t$peak_height[i])
      )
    } else {
      merged <- append(merged, list(t[i,]))
    }
  }
  
  t_merged <- do.call(rbind, merged) |> 
    mutate(peak_group = row_number()) |> 
    select(peak_group, everything())
  
  return(t_merged)
}

Green <- allLong |> 
  dplyr::filter(species == 'green', enzyme == 'CYP') |> 
  normalize_byspecies() |> 
  mutate(time_point = factor(time_point, levels = c(0,2,5,15,30))) |> 
  drop_na()

peak_df <- Green |>
  group_by(sample, time_point, species) |>
  group_modify(~ peak_ranges(.x)) |>
  ungroup()

collapsed_peak_ranges <- peak_df |>
  group_by(sample, species) |>       
  arrange(apex_time) |>
  mutate(
    peak_group = cumsum(c(TRUE, diff(apex_time) > 0.15))
  ) |>
  group_by(sample, species, peak_group) |>
  summarize(
    start_time = min(start_time),
    apex_time  = apex_time[which.max(peak_height)],
    end_time   = max(end_time),
    peak_height = max(peak_height),
    .groups = "drop"
  )

# filter for plotting just m4
green_peaks <- collapsed_peak_ranges |> 
  dplyr::filter(sample == "m4") |> 
  dplyr::select(-c(1,2))


p <- ggplot(Green, aes(time, intensity)) +
  geom_line(alpha = 0.6) +
  geom_rect(
    data = green_peaks,
    aes(xmin = start_time, xmax = end_time,
        ymin = -Inf, ymax = Inf),
    fill = "red",
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = green_peaks,
    aes(x = apex_time, y = peak_height, label = peak_group),
    inherit.aes = FALSE,
    fill = "white",
    alpha = 0.8,
    label.size = 0.2
  ) +
  labs(
    title = "Green Species"
  )+
  theme_bw()
  
p

ggsave(
  filename = "images/Peaks_Green.png",
  plot = p,
    dpi = 300,
    width = 8,
    height = 6,
    units = "in"
  )

#ggplotly(p, tooltip = "text", dynamicTicks = TRUE)