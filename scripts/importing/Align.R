ptw_function <- function(df) {
  meta <- df |> 
    dplyr::select(time_point, species, sample, enzyme, instrument)

  Num <- df |> 
    dplyr::select(where(is.numeric))

  reference <- df |> 
    dplyr::filter(time_point == 30, sample == "m4", enzyme == "CYP+UGT") |> 
    dplyr::slice(1)

  samples <- df |> 
    dplyr::filter(!(time_point == 30 & sample == "m4" & enzyme == "CYP+UGT"))

  reference_meta <- reference |> 
    dplyr::select(time_point, species, sample, enzyme, instrument)

  samples_meta <- samples |> 
    dplyr::select(time_point, species, sample, enzyme, instrument)

  reference_num <- reference |> 
    dplyr::select(where(is.numeric))

  samples_num <- samples |> 
    dplyr::select(where(is.numeric))
  
  res <- ptw(reference_num, samples_num, warp.type = "individual")

  warped_df <- as.data.frame(res$warped.sample)
  reference_df <- as.data.frame(res$reference)

  colnames(warped_df) <- colnames(reference_df)

  warped_full <- bind_cols(samples_meta, warped_df)
  reference_full <- bind_cols(reference_meta, reference_df)

  final_df <- bind_rows(reference_full, warped_full)

  return(final_df)
}

library(NoSleepR)
nosleep_on()

allWide <- correct_baseline(allWide_un, meta_in = TRUE)

species_df = split_by_species(allWide)

green <- species_df[['green']]

resGreen <- ptw_function(green)
write_parquet(resGreen, "data/processed/greenPTW_withinSpecies.parquet")

white <- species_df[['white']]

resWhite <- ptw_function(white)
write_parquet(resWhite, "data/processed/greenPTW_withinSpecies.parquet")

early <- species_df[['early']]

resEarly <- ptw_function(early)
write_parquet(resEarly, "data/processed/greenPTW_withinSpecies.parquet")

wyoming <- species_df[['wyoming']]

resWyoming <- ptw_function(wyoming)
write_parquet(resWyoming, "data/processed/greenPTW_withinSpecies.parquet")

aligned_species <- bind_rows(resWyoming, resWhite, resGreen, resEarly)
write_parquet(aligned_species, "data/processed/allPTW_withinSpecies.parquet")

Wyoming <- wide_to_long(resWyoming)
Early <- wide_to_long(resEarly)
White <- wide_to_long(resWhite)
Green <- wide_to_long(resGreen)


combined <- bind_rows(Wyoming, Early, White, Green)

write_csv(combined, "data/processed/alignedSpecies.csv")

nosleep_off()

plot_by_time(combined)
allLong <- combined
