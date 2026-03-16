library(plotly)
library(tidyverse)

combinedLong <- combinedLong |> 
  mutate(time = as.numeric(time)) |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30")))

LongPlot <- combinedLong |>
  filter(instrument == "DAD", species == 'green') |>
  mutate(sample_id = paste(enzyme, sample, time_point, species, sep = "_")) |>
  ggplot(
    aes(
      x = time,
      y = intensity,
      color = time_point,
      group = sample_id,
      text = paste(
        "Sample:", sample_id,
        "<br>Time:", round(time, 3),
        "<br>Intensity:", signif(intensity, 4),
        "<br>Species:", species
      )
    )
  ) +
  geom_line()

ggplotly(LongPlot, tooltip = "text", dynamicTicks = TRUE)


Wide <-  correct_baseline(combinedWide, meta_in = TRUE)
correctedLong <- wide_to_long(Wide) 
correctedLong <- correctedLong |> 
  mutate(time = as.numeric(time)) |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30")))


LongPlot_baseline <- correctedLong |>
  filter(instrument == "DAD", species == 'wyoming') |>
  mutate(sample_id = paste(enzyme, sample, time_point, species, sep = "_")) |>
  ggplot(
    aes(
      x = time,
      y = intensity,
      color = time_point,
      group = sample_id,
      text = paste(
        "Sample:", sample_id,
        "<br>Time:", round(time, 3),
        "<br>Intensity:", signif(intensity, 4),
        "<br>Species:", species
      )
    )
  ) +
  geom_line()

ggplotly(LongPlot_baseline, tooltip = "text", dynamicTicks = TRUE)

