library(pracma)
library(tidyverse) # for data manipulation and visualization
library(patchwork) # for combining plots
library(here) # for file path management
library(arrow)

t1 <- combinedWide |> 
  select(-c(enzyme,instrument,species,time_point,sample)) |> 
  as.matrix()


original <- t1 |> 
  as.data.frame() |> 
  pivot_longer(
    everything()) |> 
  mutate(index = row_number(),
        name = as.numeric(name))



smoothed <- smooth(combinedWide)

ggplot()+
  geom_point(data = original, aes(x = index, y = value), color = 'blue', alpha = .01, size = .1 )+
  geom_point(data = smoothed, aes(x = index, y = t2), color = 'red', alpha = .01, size = .1)
