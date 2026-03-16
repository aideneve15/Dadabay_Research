"
This script is an outline for calculating 2 fold change, p-values
And plotting with volcano style plots

The 2 fold change is between time points. 
Starting by doing the 2 fold change compared to time point 0
"
white_peak_areas_CYP <- white_peak_areas |> 
  dplyr::filter(enzyme == 'CYP')

library(dplyr)
library(broom)
library(purrr)

white_volcano_stats <- white_peak_areas_CYP |> 
  group_by(enzyme, species, peak_number, time_point) |> 
  filter(time_point != 0) |> 
  group_modify(~{
    
    # extract baseline (time_point == 0) for same enzyme/species/peak_number
    baseline <- white_peak_areas |> 
      filter(
        enzyme == .y$enzyme,
        species == .y$species,
        peak_number == .y$peak_number,
        time_point == 0
      ) |> 
      pull(area)
    
    # run t.test for this time_point
    test <- t.test(.x$area, baseline)
    
    tibble(
      time_point = unique(.x$time_point),
      mean_area = mean(.x$area),
      mean_baseline = mean(baseline),
      fc = mean(.x$area) / mean(baseline),
      p_value = test$p.value
    )
    
  }) |> 
  ungroup()

white_peak_areas_with_stats <- white_peak_areas |> 
  left_join(
    white_volcano_stats,
    by = c("enzyme", "species", "peak_number", "time_point")
  ) |> 
  mutate(fc = ifelse(time_point == 0, 1, fc))

white_peak_areas_with_stats <- white_peak_areas_with_stats |> 
  mutate(log2fc = log2(fc),
        log10p_value = -log10(p_value),
      peak_number = factor(peak_number),
      time_point = factor(time_point, levels = c(0,2,5,15,30)))

white_peak_areas_with_stats <- white_peak_areas_with_stats |>
  dplyr::mutate(rt_label = paste0(round(apex_time, 2), " RT ", time_point," min"))

white_peak_areas_with_stats$diffexpressed <- "Insignificant"
white_peak_areas_with_stats$diffexpressed[white_peak_areas_with_stats$log2fc > 0.4 & white_peak_areas_with_stats$p_value < 0.05] <- "Up Regulated"
white_peak_areas_with_stats$diffexpressed[white_peak_areas_with_stats$log2fc < -0.4 & white_peak_areas_with_stats$p_value < 0.05] <- "Down Regulated"

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("Down Regulated", "Up Regulated", "Insignificant")

white_peak_areas_with_stats$delabel <- NA
white_peak_areas_with_stats$delabel[white_peak_areas_with_stats$diffexpressed != "Insignificant"] <- white_peak_areas_with_stats$rt_label[white_peak_areas_with_stats$diffexpressed != "Insignificant"]

plot_df <- white_peak_areas_with_stats |>
  dplyr::group_by(species, enzyme, time_point, peak_number) |>
  dplyr::slice(1)

library(ggrepel)
p <- ggplot(data = plot_df, 
  aes(x = log2fc, y = log10p_value, color = diffexpressed)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.4, 0.4), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  geom_text_repel(aes(label = delabel), show.legend = FALSE) +
  theme_minimal() +
  scale_colour_manual(values = c("grey", "blue", "red")) +
  labs(
    title = "White Sagebrush with CYP Enzyme",
    subtitle = "Fold Change Compared to Zero Time Point",
    x = "log2(FoldChange)",
    y = "-log10(p-value)",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

ggsave('./images/volcano/WhiteCYP.png', p, width = 6, height = 4, units = 'in', dpi = 600)
