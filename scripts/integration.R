## Take two RT time values and integrate between. Create kinetic graphs
## Compare the kinetic graphs before and after baseline 

## Create a sample added to everything
library(tidyverse)

combinedLong <- combinedLong |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30")))

rt1 = 7.5
rt2 = 7.7

# Find rt value closest to the rt supplied
idx1 <- which.min(abs(combinedLong$time - rt1))
idx2 <- which.min(abs(combinedLong$time - rt2))

combinedLong[idx1, ]
combinedLong[idx2, ]


integrate <- combinedLong |> 
  slice(idx1:idx2) |> 
  mutate(
    dt = time - lag(time),
    intensity_mid = (intensity + lag(intensity)) / 2,
    area = dt * intensity_mid
  ) |> 
  dplyr::filter(!is.na(area))

sum(integrate$area)

ggplot(integrate) +
  geom_line(aes(x = time, y = intensity), linewidth = 1) +
  geom_rect(
    aes(
      xmin = time - dt,
      xmax = time,
      ymin = 0,
      ymax = intensity_mid
    ),
    fill = "blue",
    alpha = 0.2
  )

riemann_integrate <- function(df, rt1, rt2){
  df |>
  dplyr::filter(time >= rt1, time <= rt2) |>
  arrange(time) |>
  mutate(
    dt = time - lag(time),
    intensity_mid = (intensity + lag(intensity)) / 2,
    area = dt * intensity_mid
  ) |> 
  dplyr::filter(!is.na(area)) |> 
  summarise(area = sum(area)) |>
  pull(area)
}

riemann_integrate(combinedLong, 7.5, 7.7)


species_df = split_by_species(combinedWide)
green <- species_df[['green']]
greenLong <- wide_to_long(green)
plot_by_time(greenLong)

greenLong <- greenLong |> 
  mutate(time = as.numeric(time))

integrated_peaks <- greenLong |>
  group_by(sample, time_point) |>
  summarise(
    peak_area = riemann_integrate(cur_data(), 5.7, 6),
    .groups = "drop"
  )

integrated_peaks |> 
  dplyr::filter(sample == 'm4') |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30"))) |> 
  arrange(time_point) |> 
  ggplot()+
  geom_point(aes(x = time_point, y = peak_area, color = time_point))+
  scale_color_viridis_d(
      option = "viridis",
      name = "Incubation Time"
    )+
  labs(
    title = "Peak Abundance Over Incubation Time.",
    subtitle = "As a function of peak area."
  )

## Plot peak lists
# arguments
# - DF = Data frame which is a peak list 
# - sample = m4, m5, m6, or average (which should display the average over all samples)
# - percent = TRUE or FALSE (convert from absolute to percent decrease)