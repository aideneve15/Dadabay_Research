folder_path <- here::here("data", "raw", "combined_early_wyoming")


# List all files in the folder
files <- list.files(path = folder_path, full.names = TRUE)

all_data_long_3 <- tibble()

# iterate through files
for(i in 1:length(files)) {
  
  # Read the CSV file
  data <- read_csv(file = files[i], col_names = c("time", "intensity")) |> 
    mutate(file = files[i])
  
  all_data_long_3 <- all_data_long_3 |> 
    bind_rows(
      data
    )
}

Long3 <- all_data_long_3 |> 
  mutate(basename = basename(file)) |> 
  dplyr::select(-file) |> 
  separate_wider_delim(
    basename,
    delim = "_",
    names = c("enzyme", "instrument", "species", "sample", "time_point")
  ) |> 
  dplyr::filter(time > 0.75 & time < 10) |> 
  mutate(
  time_point = str_extract(time_point, "^[^.]+"),
  sample = tolower(sample)
)

folder_path <- here::here("data", "raw", "combined_white_green")


# List all files in the folder
files <- list.files(path = folder_path, full.names = TRUE)

all_data_long_2 <- tibble()

# iterate through files
for(i in 1:length(files)) {
  
  # Read the CSV file
  data <- read_csv(file = files[i], col_names = c("time", "intensity")) |> 
    mutate(file = files[i])
  
  all_data_long_2 <- all_data_long_2 |> 
    bind_rows(
      data
    )
}

Long2 <- all_data_long_2 |> 
  mutate(basename = basename(file)) |> 
  dplyr::select(-file) |> 
  separate_wider_delim(
    basename,
    delim = "_",
    names = c("enzyme", "instrument", "species", "sample", "time_point")
  ) |> 
  dplyr::filter(time > 0.75 & time < 10) |> 
  mutate(
  time_point = str_extract(time_point, "^[^.]+"),
  sample = tolower(sample)
)


CombinedLong <- bind_rows(Long2, Long3)

'{
all_times <- CombinedLong_in |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

CombinedLong <- CombinedLong_in |> 
  group_by(enzyme, sample, instrument, species, time_point) |>
  group_modify(~ {
    interp_results <- approx(x = .x$time, y = .x$intensity, xout = all_times$time)
    tibble(
      time = interp_results$x,
      intensity = interp_results$y
    )
  }) |> 
  ungroup()


na_times <- CombinedLong |> 
  dplyr::filter(is.na(intensity)) |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

CombinedLong <- CombinedLong |> 
  dplyr::filter(! time %in% na_times$time) |> 
  mutate(time_point = as.character(time_point))

CombinedWide <- CombinedLong |> 
  mutate(time_point = as.character(time_point)) |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  ) |> 
  dplyr::filter(instrument == "DAD")

library(arrow)

write_parquet(CombinedWide, "data/processed/CombinedWide.parquet")
write_parquet(CombinedLong, "data/processed/CombinedLong.parquet")

}'