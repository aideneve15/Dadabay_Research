"
Load the Long data frame, we use the baseline corrected and aligned.
"

allLong <- read_csv('data/processed/alignedSpecies.csv')

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
GreenCYP_peaks <- collapsed_peak_ranges 


Green <- allLong |> 
  dplyr::filter(species == 'green', enzyme == 'UGT') |> 
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
GreenUGT_peaks <- collapsed_peak_ranges 


Green <- allLong |> 
  dplyr::filter(species == 'green', enzyme == 'CYP+UGT') |> 
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

GreenCom_peaks <- collapsed_peak_ranges |> 
  dplyr::filter(start_time > 0.8) 


p <- ggplot(Green, aes(time, intensity)) +
  geom_line(alpha = 0.6) +
  geom_rect(
    data = GreenCom_peaks,
    aes(xmin = start_time, xmax = end_time,
        ymin = -Inf, ymax = Inf),
    fill = "red",
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = GreenCom_peaks,
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