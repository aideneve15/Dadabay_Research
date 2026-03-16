Long_in <- read_csv('data/processed/fullAligned.csv') 


Long <- Long_in |> 
  normalize_byspecies() |> 
  mutate(time_point = factor(time_point, levels = c(0,2,5,15,30))) 

Wide <- long_to_wide(Long) 

species_df <- split_by_species(Wide)

Green <- Long |> dplyr::filter(species == 'green')
#sample_df <- Long |> dplyr::filter(species == 'green', time_point == '0', )

# Now we have baseline corrected, aligned, and normalized samples.

detect_peaks <- function(data, deriv_threshold = 6e-1, intensity_threshold = 0.25) {
  data |> 
    group_by(sample, time_point, species) |> 
    group_modify(~ {
      
      df <- .x |> arrange(time)
      
      time <- df$time
      intensity <- df$intensity
      
      # First & second derivatives
      df$first_derivative  <- gradient(intensity, time)
      df$second_derivative <- gradient(df$first_derivative, time)
      
      d_sign <- sign(df$first_derivative)
      
      # Candidate peaks (+ to - slope change)
      peak_indices_in <- which(diff(d_sign) == -2) + 1
      
      # Apply intensity filter
      peak_indices <- peak_indices_in[
        df$intensity[peak_indices_in] > intensity_threshold
      ]
      
      if (length(peak_indices) == 0) {
        return(tibble())   # return empty if no peaks
      }
      
      peak_indices <- peak_indices_in[df$intensity[peak_indices_in] > intensity_threshold]
      peak_times <- df$time[peak_indices]
      peak_intensities <- df$intensity[peak_indices]

      peak_df <- data.frame(
        time = peak_times,
        intensity = peak_intensities,
        second_derivative = df$second_derivative[peak_indices]
      )
      # Keep only concave-down peaks
      peak_times <- peak_df |> dplyr::filter(second_derivative < 0)
      
      peak_df
    }) |> 
    ungroup()
}

peak_times <- detect_peaks(Green)

library(dplyr)

get_peak_ranges <- function(df, peak_times, threshold_fraction = 0.10) {
  
  merge_ranges <- function(ranges_df) {
    if (nrow(ranges_df) == 0) return(ranges_df)
    
    ranges_df <- ranges_df[order(ranges_df$start), ]
    merged <- list()
    current <- ranges_df[1, ]
    
    for (i in 2:nrow(ranges_df)) {
      if (ranges_df$start[i] <= current$end) {
        current$end <- max(current$end, ranges_df$end[i])
      } else {
        merged <- append(merged, list(current))
        current <- ranges_df[i, ]
      }
    }
    
    merged <- append(merged, list(current))
    do.call(rbind, merged)
  }
  
  peak_times |>
    group_by(species, time_point, sample) |>
    group_modify(~ {
      
      group_peaks <- .x
      group_key   <- .y
      
      # Filter df using group key values
      group_data <- df |>
        dplyr::filter(
          species    == group_key$species,
          time_point == group_key$time_point,
          sample     == group_key$sample
        ) |>
        arrange(time)
      
      if (nrow(group_data) == 0) return(NULL)
      
      peak_ranges <- lapply(group_peaks$time, function(peak_time) {
        
        peak_row <- which.min(abs(group_data$time - peak_time))
        peak_intensity <- group_data$intensity[peak_row]
        threshold <- peak_intensity * threshold_fraction
        
        start_row <- peak_row
        while (start_row > 1 &&
               group_data$intensity[start_row] > threshold) {
          start_row <- start_row - 1
        }
        
        end_row <- peak_row
        while (end_row < nrow(group_data) &&
               group_data$intensity[end_row] > threshold) {
          end_row <- end_row + 1
        }
        
        data.frame(
          start = group_data$time[start_row],
          apex  = group_data$time[peak_row],
          end   = group_data$time[end_row]
        )
      })
      
      peak_ranges_df <- bind_rows(peak_ranges)
      merge_ranges(peak_ranges_df)
      
    }) |>
    ungroup()
}

peak_ranges_merged <- get_peak_ranges(Green, peak_times, threshold_fraction = 0.10)


# peak_ranges <- peak_ranges_merged |>
#   group_by(time_point) |>
#   arrange(apex) |>
#   mutate(peak_id = cumsum(c(TRUE, diff(apex) > 0.2))) |>  # adjust tolerance
#   group_by(time_point, peak_id) |>
#   summarise(
#     start = min(start),
#     apex  = mean(apex),
#     end   = max(end),
#     .groups = "drop"
#   )


merge_similar_peaks <- function(df, tolerance = 0) {
  
  library(dplyr)
  
  # Sort by start time
  df <- df |> arrange(start)
  
  if (nrow(df) == 0) return(df)
  
  merged <- list()
  current <- df[1, ]
  
  for (i in 2:nrow(df)) {
    
    # Check overlap (with optional tolerance)
    if (df$start[i] <= current$end + tolerance) {
      
      # Expand range
      current$start <- min(current$start, df$start[i])
      current$end   <- max(current$end, df$end[i])
      current$apex  <- mean(c(current$apex, df$apex[i]))
      
    } else {
      
      merged[[length(merged) + 1]] <- current
      current <- df[i, ]
    }
  }
  
  merged[[length(merged) + 1]] <- current
  
  dplyr::bind_rows(merged)
}

peak_ranges_final_in <- merge_similar_peaks(peak_ranges, tolerance = 0.02)
peak_ranges_final <- peak_ranges_final_in



p <- ggplot(Long, aes(x = time, y = intensity)) +
  geom_line(aes(color = time_point)) +
  theme_minimal() +
  labs(title = "Chromatogram with Peak Ranges",
       x = "Retention Time",
       y = "Intensity")

# Add peak ranges as shaded rectangles
for (i in 1:nrow(peak_ranges_final)) {
  p <- p + 
    annotate("rect",
             xmin = peak_ranges_final$start[i],
             xmax = peak_ranges_final$end[i],
             ymin = -Inf,
             ymax = Inf,
             alpha = 0.2, fill = "red") +
    annotate("segment",
             x = peak_ranges_final$apex[i],
             xend = peak_ranges_final$apex[i],
             y = 0,
             yend = max(sample_df$intensity),
             color = "darkred", linetype = "dashed")
}

# Show plot
print(p)


results <- list()
Green <- Green |> mutate(sample_id = paste(species, sample, time_point, sep = '_'))
samples <- unique(Green$sample_id)

# Loop over samples
for (s in samples) {
  
  # Subset for this sample
  sample_data <- Green |> dplyr::filter(sample_id == s)
  
  # Loop over peak ranges
  for (i in 1:nrow(peak_ranges_final)) {
    
    peak <- peak_ranges_final[i, ]
    
    # Apply Riemann integration
    area <- riemann_integrate(
      df = sample_data,
      rt1 = peak$start,
      rt2 = peak$end
    )
    
    # Save result
    results[[length(results) + 1]] <- data.frame(
      sample = s,
      peak_id = peak$peak_id,
      area = area
    )
  }
}

# Combine all results into a dataframe
integrated_peaks <- bind_rows(results)


peak_wide <- integrated_peaks |> 
  pivot_wider(
    id_cols = sample, 
    names_from = peak_id, 
    values_from = area
  )
