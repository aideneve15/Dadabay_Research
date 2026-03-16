species_df = split_by_species(allWide)
green <- species_df[['green']]
white <- species_df[['white']]
early <- species_df[['early']]
wyoming <- species_df[['wyoming']]

green0 <- green |> 
  filter(time_point == 0, sample == 'm4', enzyme = 'CYP') |> 
  select(where(is.numeric))

white0 <- white |> 
  filter(time_point == 0, sample == 'm4', enzyme = 'CYP') |> 
  select(where(is.numeric))

early0 <- early |> 
  filter(time_point == 0, sample == 'm4', enzyme = 'CYP') |> 
  select(where(is.numeric))

wyoming0 <- wyoming |> 
  filter(time_point == 0, sample == 'm4', enzyme = 'CYP') |> 
  select(where(is.numeric))

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
    select(time_point, species, sample, enzyme, instrument)

sample <- Num |>
  select(species, sample, time_point, time, intensity) |> 
  pivot_wider(
    id_cols = c(species, sample, time_point, enzyme),
    names_from = time,
    values_from = intensity)

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

warpedLong <- warped_df |>  
  pivot_longer(
    cols = everything(),
    names_to = "time",
    values_to = "intensity")


# ------------------------------
# 3. Plot
# ------------------------------
ggplot() +
  geom_line(data = warped_df, 
            aes(x = time, y = intensity, color = species, group = species), size = 1) +
  geom_line(data = reference_dfLong, 
            aes(x = time, y = intensity), color = 'black', size = 1) +
  facet_wrap(~ species) +
  theme_minimal() +
  labs(title = "Warped spectra by species",
       x = "Time",
       y = "Intensity")


low <- warped_df |> 
  filter(is.na(intensity), time < 5) |>
  pull(time) |> 
  max()

high <- warped_df |> 
  filter(is.na(intensity), time > 5) |> 
  pull(time) |> 
  min()



warped_clean <- warped_df|> filter(time > low, time < high)
reference_clean <- reference_dfLong |> filter(time > low, time < high)

warped_list <- warped_clean |> 
  mutate(sample_id = paste0(species,time_point,sample)) |> 
  group_by(sample_id) |> 
  group_split()

library(purrr)

corsB <- map_dbl(
  warped_list,
  ~ cor(.x$intensity, 
        reference_clean$intensity)
)

corsB

# After warping there are na's, however they are only at the beginning and end 
# this justifies the dropping of na's
{
warped_df |> 
  filter(is.na(intensity)) |> 
  distinct(time)

warped_df |> 
  filter(is.na(intensity), time > 5) |> 
  summarise(mean_time = mean(time, na.rm = TRUE))

warped_df |> 
  filter(is.na(intensity), time < 5) |> 
  summarise(mean_time = mean(time, na.rm = TRUE))
}