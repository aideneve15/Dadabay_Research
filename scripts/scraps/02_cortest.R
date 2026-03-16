library(purrr)

warped_dfB <- read_csv("data/processed/fullAligned.csv")
reference_dfLongB <- read_csv("data/processed/fullAligned_reference.csv")

low <- warped_dfB |> 
  filter(is.na(intensity), time < 5) |>
  pull(time) |> 
  max()

high <- warped_dfB |> 
  filter(is.na(intensity), time > 5) |> 
  pull(time) |> 
  min()



warped_cleanB <- warped_dfB |> filter(time > low, time < high)
reference_cleanB <- reference_dfLongB |> filter(time > low, time < high)

warped_listB <- warped_cleanB |> 
  mutate(sample_id = paste0(species,time_point,sample)) |> 
  group_by(sample_id) |> 
  group_split()

corsB <- map_dbl(
  warped_listB,
  ~ cor(.x$intensity, 
        reference_cleanB$intensity)
)

corsB

warped_dfA <- read_csv("data/processed/alignedSpecies.csv")

lowA <- warped_dfA |> 
  filter(is.na(intensity), time < 5) |>
  pull(time) |> 
  max()

highA <- warped_dfA |> 
  filter(is.na(intensity), time > 5) |> 
  pull(time) |> 
  min()



warped_cleanA <- warped_dfA |> filter(time > lowA, time < highA)
reference_cleanA <- reference_dfLongB |> filter(time > lowA, time < highA)

warped_listA <- warped_cleanA |> 
  mutate(sample_id = paste0(species,time_point,sample)) |> 
  group_by(sample_id) |> 
  group_split()

corsA <- map_dbl(
  warped_listA,
  ~ cor(.x$intensity, 
        reference_cleanA$intensity)
)


mean(corsA)
var(corsA)

mean(corsB)
var(corsB)
