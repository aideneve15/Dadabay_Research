allLong_in <- bind_rows(CombinedLong, UGTLong, CYPLong) |> mutate(sampleID = paste(enzyme, species, sample, time_point, sep = "_"))

all_times <- allLong_in |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

interpolatedLong <- allLong_in |> 
  group_by(enzyme, sample, instrument, species, time_point) |>
  group_modify(~ {
    interp_results <- approx(x = .x$time, y = .x$intensity, xout = all_times$time)
    tibble(
      time = interp_results$x,
      intensity = interp_results$y
    )
  }) |> 
  ungroup()


na_times <- interpolatedLong |> 
  dplyr::filter(is.na(intensity)) |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

allLong_un <- interpolatedLong |> 
  dplyr::filter(! time %in% na_times$time) |> 
  mutate(time_point = as.character(time_point))

allWide_un <- allLong_un |> 
  mutate(time_point = as.character(time_point)) |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  ) |> 
  dplyr::filter(instrument == "DAD")

