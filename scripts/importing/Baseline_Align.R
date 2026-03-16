library(NoSleepR)
nosleep_on()

### This script was for alignment to a single sample ###

allWide <- correct_baseline(allWide_un, meta_in = TRUE)

species_df <- split(allWide, allWide[['species']])

species_df = split_by_species(allWide)
green <- species_df[['green']]
white <- species_df[['white']]
early <- species_df[['early']]
wyoming <- species_df[['wyoming']]

green0 <- green |> 
  dplyr::filter(time_point == 0, sample == 'm4', enzyme == 'CYP') |> 
  dplyr::select(where(is.numeric))

white0 <- white |> 
  dplyr::filter(time_point == 0, sample == 'm4', enzyme == 'CYP') |> 
  dplyr::select(where(is.numeric))

early0 <- early |> 
  dplyr::filter(time_point == 0, sample == 'm4', enzyme == 'CYP') |> 
  dplyr::select(where(is.numeric))

wyoming0 <- wyoming |> 
  dplyr::filter(time_point == 0, sample == 'm4', enzyme == 'CYP') |> 
  dplyr::select(where(is.numeric))

main0 = green0 + white0 + early0 + wyoming0
mainLong <- wide_to_long(main0)

mainlong <- normalize(mainLong)

reference <- mainlong |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  )

Num <- allWide |> group_by(species) |> 
  normalize_byspecies(from_wide = TRUE)

Num <- Num |> 
  mutate(
    species = as.character(species),
    sample = as.character(sample),
    time_point = as.character(time_point)
  )


meta <- Num |>
    dplyr::select(time_point, species, sample, enzyme, instrument)

sample <- Num |>
  dplyr::select(species, sample, time_point, time, intensity, enzyme) |> 
  pivot_wider(
    id_cols = c(species, sample, time_point, enzyme),
    names_from = time,
    values_from = intensity)

sample_meta <- sample |>
  dplyr::select(species, sample, time_point, enzyme)

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

colnames(warped_df) <- colnames(reference_df)
warped_df <- dplyr::bind_cols(sample_meta, warped_df)
write_parquet(warped_df, "data/processed/warped_df_singleReference.parquet")


warpedLong <- wide_to_long(warped_df) |> mutate(time = as.numeric(time))
write_parquet(warpedLong, "data/processed/warpedLong_singleReference.parquet")


reference_dfLong <- reference_df |> 
  pivot_longer(
    cols = everything(),
    names_to = 'time',
    values_to = 'intensity'
  ) |> 
  mutate(time = as.numeric(time))


# ------------------------------
# 3. Plot
# ------------------------------
ggplot() +
  geom_line(data = warpedLong, 
            aes(x = time, y = intensity, color = species, size = 1)) +
  geom_line(data = reference_dfLong, 
            aes(x = time, y = intensity), color = 'black', size = 1) +
  theme_minimal() +
  labs(title = "Warped spectra by species",
       x = "Time",
       y = "Intensity")


low <- warpedLong |> 
  dplyr::filter(is.na(intensity), time < 5) |>
  pull(time) |> 
  max()

high <- warpedLong |> 
  dplyr::filter(is.na(intensity), time > 5) |> 
  pull(time) |> 
  min()

warpedLong %>%
  group_by(species, sample, enzyme, time_point) %>%
  summarise(
    na_count = sum(is.na(intensity)),
    .groups = "drop"
  )

warped_clean <- warped_df|> dplyr::filter(time > low, time < high)
reference_clean <- reference_dfLong |> dplyr::filter(time > low, time < high)

warped_list <- warped_clean |> 
  mutate(sample_id = paste0(species,time_point,sample,enzyme)) |> 
  group_by(sample_id) |> 
  group_split()

