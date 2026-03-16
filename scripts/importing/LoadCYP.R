library(tidyverse) # for data manipulation and visualization
library(patchwork) # for combining plots
library(here) # for file path management
library(arrow) # for saving dataframes
library(ptw)
library(pracma)
library(ggforce)
library(purrr)


# Set the path to your folder
folder_path <- here::here("data", "raw", "CYP_green_white")


# List all files in the folder
files <- list.files(path = folder_path, full.names = TRUE)

all_data_long_1 <- tibble()

# iterate through files
for(i in 1:length(files)) {
  
  # Read the CSV file
  data <- read_csv(file = files[i], col_names = c("time", "intensity")) |> 
    mutate(file = files[i])
  
  all_data_long_1 <- all_data_long_1 |> 
    bind_rows(
      data
    )
}


Long1 <- all_data_long_1 |> 
  mutate(basename = basename(file)) |> 
  dplyr::select(-file) |> 
  separate_wider_delim(
    basename,
    delim = "_",
    names = c("enzyme", "instrument", "species", "time_point", "sample")
  ) |> 
  dplyr::filter(time > 0.75 & time < 10) |> 
  mutate(
  sample = str_extract(sample, "^[^.]+")
)


folder_path <- here::here("data", "raw", "CYP_early_wyoming")


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


CYPLong <- bind_rows(Long2, Long1)

# allLong |> 
#   filter(instrument == "DAD") |> 
#   group_by(species) |> 
#   ggplot(aes(x=time, y = intensity, group = 1, color = time_point))+
#   geom_line()+
#   facet_wrap(~species)

'{
all_times <- CYPLong_in |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

CYPLong <- CYPLong_in |> 
  group_by(enzyme, sample, instrument, species, time_point) |>
  group_modify(~ {
    interp_results <- approx(x = .x$time, y = .x$intensity, xout = all_times$time)
    tibble(
      time = interp_results$x,
      intensity = interp_results$y
    )
  }) |> 
  ungroup()


na_times <- CYPLong |> 
  dplyr::filter(is.na(intensity)) |> 
  dplyr::select(time) |> 
  distinct() |> 
  arrange(time)

CYPLong <- CYPLong |> 
  dplyr::filter(! time %in% na_times$time) |> 
  mutate(time_point = as.character(time_point))

# bin_width <- 0.001

# binned_data <- allLong |>
#   mutate(time_bin = floor(time / bin_width) * bin_width) |>
#   group_by(enzyme, sample, instrument, species, time_point, time_bin) |>
#   summarise(
#     intensity = sum(intensity, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   rename(time = time_bin)

CYPWide <- CYPLong |> 
  mutate(time_point = as.character(time_point)) |> 
  pivot_wider(
    names_from = time,
    values_from = intensity
  ) |> 
  dplyr::filter(instrument == "DAD")

library(arrow)

write_parquet(CYPWide, "data/processed/CYPWide.parquet")
write_parquet(CYPLong, "data/processed/CYPLong.parquet")

}'