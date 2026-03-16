#allLong <- warped_df
allLong <- wide_to_long(allWide) |> mutate(time = as.numeric(time))

library(arrow)
write_parquet(allWide, "data/processed/allLong_processed.parquet")
write_parquet(allLong, "data/processed/allWide_processed.parquet")
