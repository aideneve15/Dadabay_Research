# peaks <- pracma::findpeaks(sample$intensity)
# peaks_index <- peaks[,2]

# peak_times <- sample$time[peaks_index]

# peaks_candidate <- sample[
#   abs(sample$first_derivative) < deriv_threshold & 
#   sample$intensity > intensity_threshold &
#   sample$second_derivative < 0, 
# ]

Long_in <- read_csv('data/processed/fullAligned.csv')
#plot_by_time(Long)

library(plotly)
library(signal)

"The Long Dataframe has been baseline corrected and aligned"
"Now normalize"

Long <- Long_in |> 
  normalize_byspecies() 

# Smoothed <- Long |> long_to_wide() |> 
#   group_by(species, time_point, sample) |> 
#   mutate(across(where(is.numeric), ~ calc_ewma(.x, a = 0.3))) |> 
#   ungroup()

# Long_df <- Smoothed |> wide_to_long()
# Long_df <- Long_df |> mutate(time = as.numeric(time))

# Long_df |>
#   group_by(species, sample, time_point) |>
#   summarise(area = riemann_integrate(cur_data(), 0.7, 10))

sample_df <- dplyr::filter(Long, sample == 'm4', species == 'early', time_point == 0)

sample_df$intensity <- sgolayfilt(sample_df$intensity, p = 3, n = 7)

time <- sample_df$time
intensity <- sample_df$intensity

# First derivative
d_intensity <- gradient(intensity, time)
sample_df$first_derivative <- d_intensity
sample_df$second_derivative <- gradient(sample_df$first_derivative, time)

deriv_threshold <- 6e-1
intensity_threshold <- 0.25

d_sign <- sign(sample_df$first_derivative)

# Candidate peaks where derivative switches from + to -
peak_indices_in <- which(diff(d_sign) == -2) + 1

# Filter by intensity threshold
peak_indices <- peak_indices_in[sample_df$intensity[peak_indices_in] > intensity_threshold]
peak_times <- sample_df$time[peak_indices]
peak_intensities <- sample_df$intensity[peak_indices]

# Optional: create a small data frame for easy inspection
peak_df <- data.frame(
  time = peak_times,
  intensity = peak_intensities,
  second_derivative = sample_df$second_derivative[peak_indices]
)

peak_times <- peak_df |> dplyr::filter(second_derivative < 0)

plot <- sample_df |> 
  dplyr::filter(instrument == "DAD") |>
  mutate(sample_id = paste(enzyme, sample, time_point, species, sep = "_")) |>
  ggplot(
    aes(
      x = time,
      y = intensity,
      color = time_point,
      group = sample_id,
  )
) +
geom_line() +
geom_vline(data = peak_times, aes(xintercept = time), size = 0.1, color = 'red') +
labs(
  title = "title_input",
  x = "Retention Time",
  y = "Intensity (au)"
)+
theme(
  plot.title = element_text(face = "bold"), 
  axis.title.x = element_text(face = "bold"),         
  axis.title.y = element_text(face = "bold"))

ggplotly(plot, tooltip = "text", dynamicTicks = TRUE)

