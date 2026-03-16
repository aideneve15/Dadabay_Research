"
This script is an outline for calculating 2 fold change, p-values
And plotting with volcano style plots

The 2 fold change is between time points. 
Starting by doing the 2 fold change compared to time point 0
"
early_peak_areas_CYP <- early_peak_areas |> 
  dplyr::filter(enzyme == 'CYP')

library(dplyr)
library(broom)
library(purrr)

early_volcano_stats <- early_peak_areas_CYP |> 
  group_by(enzyme, species, peak_number, time_point) |> 
  filter(time_point != 0) |> 
  group_modify(~{
    
    # extract baseline (time_point == 0) for same enzyme/species/peak_number
    baseline <- early_peak_areas |> 
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

early_peak_areas_with_stats <- early_peak_areas |> 
  left_join(
    early_volcano_stats,
    by = c("enzyme", "species", "peak_number", "time_point")
  ) |> 
  mutate(fc = ifelse(time_point == 0, 1, fc))

early_peak_areas_with_stats <- early_peak_areas_with_stats |> 
  mutate(log2fc = log2(fc),
        log10p_value = -log10(p_value),
      peak_number = factor(peak_number),
      time_point = factor(time_point, levels = c(0,2,5,15,30)))

early_peak_areas_with_stats <- early_peak_areas_with_stats |>
  dplyr::mutate(rt_label = paste0(round(apex_time, 2), " RT ", time_point," min"))

early_peak_areas_with_stats$diffexpressed <- "Insignificant"
early_peak_areas_with_stats$diffexpressed[early_peak_areas_with_stats$log2fc > 0.4 & early_peak_areas_with_stats$p_value < 0.05] <- "Up Regulated"
early_peak_areas_with_stats$diffexpressed[early_peak_areas_with_stats$log2fc < -0.4 & early_peak_areas_with_stats$p_value < 0.05] <- "Down Regulated"

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("Down Regulated", "Up Regulated", "Insignificant")

early_peak_areas_with_stats$delabel <- NA
early_peak_areas_with_stats$delabel[early_peak_areas_with_stats$diffexpressed != "Insignificant"] <- early_peak_areas_with_stats$rt_label[early_peak_areas_with_stats$diffexpressed != "Insignificant"]

plot_df <- early_peak_areas_with_stats |>
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
  scale_colour_manual(values = c("blue", "grey", "red")) +
  labs(
    title = "Early Sagebrush with UGT Enzyme",
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

ggsave('./images/volcano/EarlyCYP.png', p, width = 6, height = 4, units = 'in', dpi = 600)
# {
# library(broom)

# early_peak_areas <- early_peak_areas |> 
#   mutate(time_point = as.numeric(as.character(time_point)))

# # Identify the last time point per enzyme/species/peak_number group
# last_timepoints <- early_peak_areas |> 
#   group_by(enzyme, species, peak_number) |> 
#   summarize(last_tp = max(time_point), .groups = "drop")

# # Calculate fold change and p-values compared to last time point
# early_volcano_stats_last <- early_peak_areas |> 
#   group_by(enzyme, species, peak_number, time_point) |> 
#   filter(time_point != last_timepoints$last_tp[match(paste(enzyme, species, peak_number),
#                                                      paste(last_timepoints$enzyme, 
#                                                            last_timepoints$species, 
#                                                            last_timepoints$peak_number))]) |> 
#   group_modify(~{
    
#     # Extract last time point value for same enzyme/species/peak_number
#     last_tp_val <- early_peak_areas |> 
#       filter(
#         enzyme == .y$enzyme,
#         species == .y$species,
#         peak_number == .y$peak_number,
#         time_point == max(early_peak_areas$time_point[early_peak_areas$enzyme == .y$enzyme &
#                                                       early_peak_areas$species == .y$species &
#                                                       early_peak_areas$peak_number == .y$peak_number])
#       ) |> 
#       pull(area)
    
#     # run t.test for this time_point vs last timepoint
#     test <- t.test(.x$area, last_tp_val)
    
#     tibble(
#       time_point = unique(.x$time_point),
#       mean_area = mean(.x$area),
#       mean_last = mean(last_tp_val),
#       fc = mean(.x$area) / mean(last_tp_val),
#       p_value = test$p.value
#     )
    
#   }) |> 
#   ungroup()

# # Merge back with original data
# early_peak_areas_with_stats_last <- early_peak_areas |> 
#   left_join(
#     early_volcano_stats_last,
#     by = c("enzyme", "species", "peak_number", "time_point")
#   ) |> 
#   mutate(fc = ifelse(time_point == max(time_point), 1, fc))  # fc = 1 for last time point

# # Add log transforms for volcano plot
# early_peak_areas_with_stats_last <- early_peak_areas_with_stats_last |> 
#   mutate(
#     log2fc = log2(fc),
#     log10p_value = -log10(p_value),
#     peak_number = factor(peak_number)
#   )

# # Volcano plot
# ggplot(data = early_peak_areas_with_stats_last, aes(x = log2fc, y = log10p_value, shape = time_point, color = peak_number)) +
#   geom_point() +
#   geom_vline(xintercept = c(-0.4, 0.4), col = "red") +
#   geom_hline(yintercept = -log10(0.05), col = "red") +
#   theme_minimal() +
#   labs(
#     title = "Early Sagebrush with CYP Enzyme (Compared to Last Timepoint)"
#   )


# # de$delabel <- NA
# # de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]


# # ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
# #     geom_point() + 
# #     theme_minimal() +
# #     geom_text()
# }