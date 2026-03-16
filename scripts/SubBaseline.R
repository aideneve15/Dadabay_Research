# Create function to section off dataframe based on RT


# Create function to get baseline 
species_df = combinedWide |> split_by_species()
green <- species_df[['green']]
green_cor <- correct_baseline(green, l = 1*10^9, maxit = 25)

greenNum <- green |> 
  select(where(is.numeric))

green_corNum <- green_cor |> 
  select(where(is.numeric))

baselineWide <- greenNum - green_corNum

baseline <- baselineWide[2,] |> 
  pivot_longer(
    cols = everything(),
    names_to = 'time',
    values_to = 'intensity'
  )

greenLong <- green[2,] |> 
  pivot_longer(
    cols = -c(species, time_point, enzyme, instrument, sample),
    names_to = 'time',
    values_to = 'intensity'
  ) |> 
  mutate(time = as.numeric(time))

baseline <- baseline |> 
  mutate(time = as.numeric(time))

ggplot() +
  geom_line(
    data = baseline,
    aes(x = time, y = intensity, color = "Baseline correction"),
    linewidth = 1
  ) +
  geom_line(
    data = greenLong,
    aes(x = time, y = intensity, color = "Raw spectra"),
    linewidth = 0.8
  ) +
  scale_color_manual(
    name = NULL,  # removes legend title (optional)
    values = c(
      "Baseline correction" = "pink",
      "Raw spectra" = "lightblue"
    )
  ) +
  labs(
    x = "Time",
    y = "Intensity"
  ) +
  theme_minimal()

### Plot iterations of the baseline
# Generate green shades (dark → light)
green_colors <- colorRampPalette(
  c("#0b3d0b", "#66c266")
)(10)

baseline_all <- map_dfr(
  1:10,
  function(i) {

    green_cor <- correct_baseline(
      green,
      l = 1e9,
      maxit = i
    )

    green_corNum <- green_cor |> 
      select(where(is.numeric))

    baselineWide <- greenNum - green_corNum

    baselineWide[2,] |> 
      pivot_longer(
        cols = everything(),
        names_to = "time",
        values_to = "intensity"
      ) |> 
      mutate(
        time = as.numeric(time),
        iterations = factor(i)
      )
  }
)

greenLong <- green[2,] |> 
  pivot_longer(
    cols = -c(species, time_point, enzyme, instrument, sample),
    names_to = "time",
    values_to = "intensity"
  ) |> 
  mutate(time = as.numeric(time))

ggplot() +
  geom_line(
    data = baseline_all,
    aes(x = time, y = intensity, color = iterations),
    linewidth = 0.5,
    alpha = 1
  ) +
  # geom_line(
  #   data = greenLong,
  #   aes(x = time, y = intensity),
  #   color = "red",
  #   alpha = 0.5,
  #   linewidth = 0.8
  # ) +
  scale_color_manual(
    values = green_colors,
    name = "Number of iterations",
    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  ) +
  labs(
    x = "Time",
    y = "Intensity"
  ) +
  theme_minimal()

ggsave("images/iterations.png")
