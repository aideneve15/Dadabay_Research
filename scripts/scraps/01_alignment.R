combinedWide <- read_parquet("data/processed/combinedWide.parquet")
combinedWide <- combinedWide |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30")))

combinedLong <- read_parquet("data/processed/combinedLong.parquet")
combinedLong <- combinedLong |> 
  mutate(time_point = factor(time_point, levels = c("0", "2", "5", "15", "30")))

species_df = split_by_species(combinedWide)
green <- species_df[['green']]
white <- species_df[['white']]
early <- species_df[['early']]
wyoming <- species_df[['wyoming']]

green0 <- green |> 
  filter(time_point == 0, sample == 'm4') |> 
  select(where(is.numeric))

white0 <- white |> 
  filter(time_point == 0, sample == 'm4') |> 
  select(where(is.numeric))

early0 <- early |> 
  filter(time_point == 0, sample == 'm4') |> 
  select(where(is.numeric))

wyoming0 <- wyoming |> 
  filter(time_point == 0, sample == 'm4') |> 
  select(where(is.numeric))

main0 = green0 + white0 + early0 + wyoming0
mainLong <- wide_to_long(main0)

mainlong <- normalize(mainLong)

reference <- mainlong |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  )

meta <- combinedWide |> group_by(species) |> 
    filter(sample == 'm4', time_point == 0) |>
    select(time_point, species, sample, enzyme, instrument)

Num <- combinedWide |> group_by(species) |> 
  filter(sample == 'm4', time_point == 0) |> 
  normalize_byspecies(from_wide = TRUE)

Num <- Num |> 
  mutate(
    species = as.character(species),
    sample = as.character(sample),
    time_point = as.character(time_point)
  )

sample <- Num |>
  select(species, sample, time_point, time, intensity) |> 
  pivot_wider(
    id_cols = c(species, sample, time_point),
    names_from = time,
    values_from = intensity) |> 
  correct_baseline(meta_in = FALSE)

reference <- reference |> correct_baseline(meta_in = FALSE)

n_ref <- ncol(reference)
n_sam <- ncol(sample)

if (n_ref > n_sam) {
  reference <- reference[, seq_len(n_sam), drop = FALSE]
} else if (n_sam > n_ref) {
  sample <- sample[, seq_len(n_ref), drop = FALSE]
}

# ------------------------------
# 1. Run PTW
# ------------------------------
res_PTWSingle <- ptw(reference, sample, warp.type = "individual")

# Extract warped samples and reference as dataframes
warped_df <- as.data.frame(res_PTWSingle$warped.sample)
reference_df <- as.data.frame(res_PTWSingle$reference)


# Make sure column names match
colnames(warped_df) <- colnames(reference_df)

# ------------------------------
# 2. Convert warped data to long format
# ------------------------------
# Number of time points in your warped spectra
n_timepoints <- ncol(warped_df)

# Species order — must match the order of rows in your 'sample' passed to PTW
species_order <- unique(Num$species)

# Create long dataframe with correct species labels
warped_dfLong <- warped_df |>
  mutate(species = species_order) |>       # assign species to each row
  pivot_longer(
    cols = -species,
    names_to = "time",
    values_to = "intensity"
  ) |>
  mutate(time = as.numeric(time),
         sample_id = paste(species, "m4", sep = "_")) 


# Similarly convert reference to long for plotting
reference_dfLong <- reference_df |>
  pivot_longer(
    cols = everything(),
    names_to = "time",
    values_to = "intensity"
  ) |>
  mutate(time = as.numeric(time))
  
write_csv(warped_dfLong, "data/processed/ptw_four_to_main_corrected.csv")
write_csv(reference_dfLong, "data/processed/ptw_four_to_main_corrected_reference.csv")
# ------------------------------
# 3. Plot
# ------------------------------
ggplot() +
  geom_line(data = warped_dfLong, 
            aes(x = time, y = intensity, color = species, group = species), size = 1) +
  geom_line(data = reference_dfLong, 
            aes(x = time, y = intensity), color = 'black', size = 1) +
  facet_wrap(~ species) +
  theme_minimal() +
  labs(title = "Warped spectra by species",
       x = "Time",
       y = "Intensity")

