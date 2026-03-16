##########--------------Coefficient of Variation------------############
"""
Signal to Noise Ratio can be defined as the reciprocal of the coefficient of variation, also its square.
To calculate the reciprocal mean/std

Should try this on baseline corrected spectra as well
"""
CofV <- combinedLong |> 
  mutate(sample_id = paste0(species, sample, time_point)) |> 
  group_by(time) |> 
  mutate(mu = mean(intensity),
        var = std(intensity)) |> 
  ungroup()

CofV |> 
  mutate(snr = var/mu) |> 
  filter(snr > 3) |>
  reframe(time, snr, sample_id) |> 
  distinct(sample_id) |> 
  arrange(sample_id) |> 
  print(n = 60)

dim(CofV)
# 1404150 Rows

#11592 Rows after filtering 

combinedFlat <- combinedWide |> correct_baseline()

flat <- wide_to_long(combinedFlat)

plot_by_time(flat)

snr_3_times <- CofV |> 
  mutate(snr = var/mu) |> 
  filter(snr > 3) |>
  reframe(time, snr, sample_id) |> 
  distinct(time)


Long <- flat |> 
    filter(instrument == "DAD") |>
    mutate(time = as.numeric(time),
      sample_id = paste(enzyme, sample, time_point, species, sep = "_")) |>
    ggplot(
      aes(
        x = time,
        y = intensity,
        color = time_point,
        group = sample_id,
    )
  ) +
  geom_line() +
  geom_vline(data = snr_3_times, aes(xintercept = time), size = 0.1, color = 'red') +
  labs(
    title = "Title",
    x = "Retention Time",
    y = "Intensity (au)"
  )+
  theme(
    plot.title = element_text(face = "bold"), 
    axis.title.x = element_text(face = "bold"),         
    axis.title.y = element_text(face = "bold"))

ggplotly(Long, tooltip = "text", dynamicTicks = TRUE)


# 1. Calculate SNR per species/time
snr_by_group <- combinedLong |>
  group_by(time, species) |>
  summarize(
    mu  = mean(intensity, na.rm = TRUE),
    std = sd(intensity, na.rm = TRUE),
    snr = mu / (std + 1e-6),
    .groups = "drop"
  )

# 2. Find time points where SNR is low across ALL species
# We look for time points where the MAX snr is still < 3
low_signal_times <- snr_by_group |>
  group_by(time) |>
  summarize(max_snr = max(snr)) |>
  filter(max_snr < 5) |>
  pull(time)

# 3. Filter the original data
combined_filtered <- combinedLong |>
  filter(!(time %in% low_signal_times))
