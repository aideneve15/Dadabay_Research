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

reference1 <- mainlong |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  )


plot(mainlong)
#Create one dataframe with the first as the combined

meta <- combinedWide |> group_by(species) |> 
    filter(sample == 'm4', time_point == 0) |>
    select(time_point, species, sample, enzyme, instrument)

Num <- combinedWide |> group_by(species) |> 
  filter(sample == 'm4', time_point == 0) |> 
  normalize(from_wide = TRUE)

sample <- Num |>
  select(species, time, intensity) |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  ) |> 
  select(-last_col()) |> 
  correct_baseline(meta_in = FALSE)

reference <- reference1 |> correct_baseline(meta_in = FALSE)

res_PTWSingle <- ptw(reference, sample, warp.type = "individual")

warped_df <- as.data.frame(res_PTWSingle$warped.sample)
reference_df <- as.data.frame(res_PTWSingle$reference)


colnames(warped_df) <- colnames(reference_df)

warped_df <- bind_cols(meta, warped_df)
write_parquet(warped_df, "data/processed/ptw_four_to_main")
write_parquet(reference_df, "data/processed/ptw_four_to_main_reference")


warped_dfLong <- wide_to_long(warped_df) |> mutate(time = as.numeric(time))
reference_dfLong <- wide_to_long(reference_df) |> mutate(time = as.numeric(time))
full <- bind_rows(warped_dfLong, reference_dfLong)

Long_df <- warped_dfLong |>  
  filter(instrument == "DAD") |>
  mutate(sample_id = paste(enzyme, sample, time_point, species, sep = "_"))

ggplot() +
  geom_line(data = Long_df, 
    aes(
        x = time,
        y = intensity,
        color = species,
        group = species,
    )
  ) +
  facet_wrap(~ species) + 
  geom_line(data = reference_dfLong, 
    aes(
        x = time,
        y = intensity
    ), color = 'black'
  )



result <- pca(warped_dfLong)
