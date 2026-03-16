# Alignment takes a significant amount of processing time 
# I will be writing into parquet

combinedWide

species_df = split_by_species(combinedWide)

green <- species_df[['green']]

greenWide <- correct_baseline(green, l = 1*10^9, meta_in = TRUE)
resGreen <- ptw_function(greenWide)

# result <- pca(resGreen, scale_na = TRUE)

# scores <- as.data.frame(result$scores)

# plot_scores(scores, PC2, PC3, title_in = "Baseline Correced Score Plot")

# write_parquet(resGreen, "data/processed/greenPTW.parquet")

# long <- wide_to_long(greenPTW)
# plot_by_time(long)

white <- species_df[['white']]

whiteWide <- correct_baseline(white, l = 1*10^9, meta_in = TRUE)
resWhite <- ptw_function(whiteWide)

write_parquet(resWhite, "data/processed/whitePTW.parquet")


early <- species_df[['early']]

earlyWide <- correct_baseline(early, l = 1*10^9, meta_in = TRUE)
resEarly <- ptw_function(earlyWide)

write_parquet(resEarly, "data/processed/earlyPTW.parquet")

wyoming <- species_df[['wyoming']]

wyomingWide <- correct_baseline(wyoming, l = 1*10^9, meta_in = TRUE)
resWyoming <- ptw_function(wyomingWide)

write_parquet(resWyoming, "data/processed/wyomingPTW.parquet")

Wyoming <- wide_to_long(resWyoming)
Early <- wide_to_long(resEarly)
White <- wide_to_long(resWhite)
Green <- wide_to_long(resGreen)


combined <- bind_rows(Wyoming, Early, White, Green)

write_csv(combined, "data/processed/alignedSpecies.csv")
