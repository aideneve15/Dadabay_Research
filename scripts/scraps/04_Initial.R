#Setup Import DataFrames

library(tidyverse) # for data manipulation and visualization
library(patchwork) # for combining plots
library(here) # for file path management
library(arrow)

dataWide <- read_parquet('data_main/processed/all_data_wide.parquet')
dataLong <- read_parquet('data_main/processed/data_long_interp.parquet')


all_data_ptw <- dataLong |> 
  select(-instrument, -enzyme, -time_point, -sample, -sb_color, -file) 


#PTW (need to split into sample and unique type of sagebrush first)
{
#Select unique files, set a seed for random but repeatable generation
#Then take the random sample from unqiue files 
files_u <- unique(dataLong$file)
set.seed(42)
files_shuffled <- sample(files_u)

mid <- ceiling(length(files_shuffled)/2)

half1_files <- files_shuffled[1:mid]
half2_files <- files_shuffled[(mid+1):length(files_shuffled)]

half1 <- dataLong |> filter(file %in% half1_files)
half2 <- dataLong |> filter(file %in% half2_files)

half1_mat <- half1 %>%
  select(file, time, intensity) %>%
  pivot_wider(names_from = time, values_from = intensity) %>%
  column_to_rownames("file") %>%   # file becomes row names
  as.matrix()

half2_mat <- half2 %>%
  select(file, time, intensity) %>%
  pivot_wider(names_from = time, values_from = intensity) %>%
  column_to_rownames("file") %>%
  as.matrix()

res <- ptw(half2_mat, half1_mat, warp.type = "global")



# Convert reference spectrum to long
ref_long <- as.data.frame(res$ref) %>%
  mutate(file = rownames(.)) %>%
  pivot_longer(-file, names_to = "time", values_to = "intensity") %>%
  mutate(time = as.numeric(time), type = "Reference")

# Convert warped spectra to long
sample_long <- as.data.frame(res$sample) %>%
  mutate(file = rownames(.)) %>%
  pivot_longer(-file, names_to = "time", values_to = "intensity") %>%
  mutate(time = as.numeric(time), type = "Warped")

# Combine
plot_data <- bind_rows(ref_long, sample_long)

ggplot(plot_data, aes(x = time, y = intensity, color = type)) +
  geom_point() +
  labs(title = "PTW Warping: Reference vs Warped Spectra") +
  theme_minimal()+
  ggsave('Dadabay_Research/images/PTW_plot.png', plot = get_last_plot, dpi = 800, create.dir = TRUE)

}
##Pre-Processing
#==================================================================================
#install.packages("ptw", type = "binary")
#library(ptw)
  #DAD
baseline_dad <- all_data_wide |> 
  filter(instrument == "DAD") |> 
  select(-instrument) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  )

metadata <- baseline_dad |> select(1:4)

# Transform the rest
transformed <- baseline_dad |> 
  select(-1, -2, -3, -4) |>   # drop metadata columns
  baseline.corr()     # apply function

# Recombine
corrected_baseline_dad <- bind_cols(metadata, transformed) 

#Uncorrected Plot
baseline_dad_plot <- baseline_dad |>                   
  pivot_longer(
    cols = -c(time_point, sb_color, sample, enzyme),  
    names_to = "rt",
    values_to = "intensity"
  ) |> 
     mutate(rt = as.numeric(rt),
            sample = as.character(sample)) |>   # convert RT column names to numeric
  ggplot(aes(x = rt, y = intensity, color = time_point)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~sample)+
  labs(
    x = "Retention Time",
    y = "Intensity",
    color = "Sample"
  ) +
  theme_minimal()

#Corrected Plot
corrected_baseline_dad_plot <- corrected_baseline_dad |>
  pivot_longer(
    cols = -c(time_point, sample, enzyme, sb_color),
    names_to = "rt",
    values_to = "intensity"
  ) |> 
     mutate(rt = as.numeric(rt),
            sample = as.character(sample)) |>   # convert RT column names to numeric
  ggplot(aes(x = rt, y = intensity, color = time_point)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~sample)+
  labs(
    x = "Retention Time",
    y = "Intensity",
    color = "Sample"
  ) +
  theme_minimal()

#MS
baseline_ms <- all_data_wide |> 
  filter(instrument == "MS") |> 
  select(-instrument) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  )

metadata_ms <- baseline_ms |> select(1:4)

# Transform the rest
transformed_ms <- baseline_ms |> 
  select(-1, -2, -3, -4) |>   # drop metadata columns
  baseline.corr()     # apply function

# Recombine
corrected_baseline_ms <- bind_cols(metadata_ms, transformed_ms) 

#Uncorrected Plot
baseline_ms_plot <- baseline_ms |>                   
  pivot_longer(
    cols = -c(time_point, sb_color, sample, enzyme),  
    names_to = "rt",
    values_to = "intensity"
  ) |> 
     mutate(rt = as.numeric(rt),
            sample = as.character(sample)) |>   # convert RT column names to numeric
  ggplot(aes(x = rt, y = intensity, color = time_point)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~sample)+
  labs(
    x = "Retention Time",
    y = "Intensity",
    color = "Incubation time"
  ) +
  theme_minimal()

#Corrected Plot
corrected_baseline_ms_plot <- corrected_baseline_ms |>
  pivot_longer(
    cols = -c(time_point, sample, enzyme, sb_color),
    names_to = "rt",
    values_to = "intensity"
  ) |> 
     mutate(rt = as.numeric(rt),
            sample = as.character(sample)) |>   # convert RT column names to numeric
  ggplot(aes(x = rt, y = intensity, color = time_point)) +
  geom_point(size = .1, alpha = .5) +
  facet_wrap(~sample)+
  labs(
    x = "Retention Time",
    y = "Intensity",
    color = "Incubation time"
  ) +
  theme_minimal()

baseline_dad_plot
corrected_baseline_dad_plot

baseline_ms_plot
corrected_baseline_ms_plot 


#PCA
{
#==================================================================================
dad_pca <- all_data_wide |> 
  filter(instrument == "DAD") |> 
  select(-instrument, -enzyme, -sb_color, -time_point, -sample) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  ) |> 
  prcomp(center = TRUE, scale. = TRUE)
  

ms_pca <- all_data_wide |> 
  filter(instrument == "MS") |> 
  select(-instrument, -enzyme, -sb_color, -time_point, -sample) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  ) |> 
  prcomp(center = TRUE, scale. = TRUE)

dad_variance <- summary(dad_pca)$importance[2, 1:15]
ms_variance <- summary(ms_pca)$importance[2, 1:15]

p1 <- tibble(PC = 1:15, Variance = dad_variance) |>
  ggplot(aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  labs(title = "Scree Plot - DAD", x = "PC", y = "Proportion of Variance")

p2 <- tibble(PC = 1:15, Variance = ms_variance) |>
  ggplot(aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  labs(title = "Scree Plot - MS", x = "PC", y = "Proportion of Variance")

}
dad_pca_baseline <- corrected_baseline_dad |> 
  select(-sb_color, -enzyme, -time_point, -sample) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  ) |> 
  prcomp(center = TRUE, scale. = TRUE)

dad_baseline_variance <- summary(dad_pca_baseline)$importance[2, 1:15]

p1_baseline <- tibble(PC = 1:15, Variance = dad_baseline_variance) |>
  ggplot(aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  labs(title = "Scree Plot - DAD, Corrected Baseline", x = "PC", y = "Proportion of Variance")


ms_pca_baseline <- corrected_baseline_ms |> 
  select(-sb_color, -enzyme, -time_point, -sample) |> 
  mutate(
    across(
      where(is.list),
      ~ purrr::map_dbl(.x, mean, na.rm = TRUE)
    )
  ) |> 
  prcomp(center = TRUE, scale. = TRUE)

ms_baseline_variance <- summary(ms_pca_baseline)$importance[2, 1:15]

p2_baseline <- tibble(PC = 1:15, Variance = ms_baseline_variance) |>
  ggplot(aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  labs(title = "Scree Plot - MS, Corrected Baseline", x = "PC", y = "Proportion of Variance")

p1 + p2
p1_baseline + p2_baseline



dad_scores <- as_tibble(dad_pca$x) %>%
  mutate(
    incubation = all_data_wide |> filter(instrument == "DAD") |> pull(time_point),
    sagebrush = all_data_wide |> filter(instrument == "DAD") |> pull(sb_color),
    enzyme = all_data_wide |> filter(instrument == "DAD") |> pull(enzyme),
    sample = all_data_wide |> filter(instrument == "DAD") |> pull(sample),
  )

p1 <- ggplot(dad_scores, aes(x = PC1, y = PC2, color = incubation, shape = sagebrush)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Biplot - DAD Measurement",
    x = "PC1",
    y = "PC2",
    color = "Incubation Time",
    shape = "Sagebrush Type"
  ) +
  theme(legend.position = "none")



dad_scores_baseline <- as_tibble(dad_pca_baseline$x) %>%
  mutate(
    incubation = all_data_wide |> filter(instrument == "DAD") |> pull(time_point),
    sagebrush = all_data_wide |> filter(instrument == "DAD") |> pull(sb_color),
    enzyme = all_data_wide |> filter(instrument == "DAD") |> pull(enzyme),
    sample = all_data_wide |> filter(instrument == "DAD") |> pull(sample),
  )

p2 <- ggplot(dad_scores_baseline, aes(x = PC1, y = PC2, color = incubation, shape = sagebrush)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Biplot - DAD Measurement after Baseline correction",
    x = "PC1",
    y = "PC2",
    color = "Incubation Time",
    shape = "Sagebrush Type"
  ) +
  theme(legend.position = "none")


ms_scores <- as_tibble(ms_pca$x) %>%
  mutate(
    incubation = all_data_wide |> filter(instrument == "DAD") |> pull(time_point),
    sagebrush = all_data_wide |> filter(instrument == "DAD") |> pull(sb_color),
    enzyme = all_data_wide |> filter(instrument == "DAD") |> pull(enzyme),
    sample = all_data_wide |> filter(instrument == "DAD") |> pull(sample),
  )


p1_ms <- ggplot(ms_scores, aes(x = PC1, y = PC2, color = incubation, shape = sagebrush)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Biplot - MS Measurement",
    x = "PC1",
    y = "PC2",
    color = "Incubation Time",
    shape = "Sagebrush Type"
  ) +
  theme(legend.position = "none")



ms_scores_baseline <- as_tibble(ms_pca_baseline$x) %>%
  mutate(
    incubation = all_data_wide |> filter(instrument == "DAD") |> pull(time_point),
    sagebrush = all_data_wide |> filter(instrument == "DAD") |> pull(sb_color),
    enzyme = all_data_wide |> filter(instrument == "DAD") |> pull(enzyme),
    sample = all_data_wide |> filter(instrument == "DAD") |> pull(sample),
  )

p2_ms <- ggplot(ms_scores_baseline, aes(x = PC1, y = PC2, color = incubation, shape = sagebrush)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Biplot - MS Measurement after Baseline correction",
    x = "PC1",
    y = "PC2",
    color = "Incubation Time",
    shape = "Sagebrush Type"
  ) +
  theme(legend.position = "none")


p1 + p2

p1_ms + p2_ms