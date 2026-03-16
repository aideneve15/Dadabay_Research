"
Peak dataframes have the following information;
peak number, start time, end time as well as an apex_time

This script should take the start and end time for all peaks and 
integrate from the corresponding species.

Group by peak sample id and peak number and integrate
"

allLong;

GreenCYP <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP', species == 'green')

GreenUGT <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'UGT', species == 'green')

GreenCom <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP+UGT', species == 'green')

GreenCom_peaks;
GreenCYP_peaks;
GreenUGT_peaks;

library(purrr)
greenCYP_peak_areas <- GreenCYP |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(GreenCYP_peaks), function(i){
      
      peak <- GreenCYP_peaks[i,]
      
      area <- riemann_integrate(
        .x,
        rt1 = peak$start_time,
        rt2 = peak$end_time
      )
      
      tibble(
        peak_number = peak$peak_group,
        start_time = peak$start_time,
        end_time = peak$end_time,
        apex_time = peak$apex_time,
        area = area
      )
      
    })
    
  }) |> 
  ungroup()

greenUGT_peak_areas <- GreenUGT |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(GreenUGT_peaks), function(i){
      
      peak <- GreenUGT_peaks[i,]
      
      area <- riemann_integrate(
        .x,
        rt1 = peak$start_time,
        rt2 = peak$end_time
      )
      
      tibble(
        peak_number = peak$peak_group,
        start_time = peak$start_time,
        end_time = peak$end_time,
        apex_time = peak$apex_time,
        area = area
      )
      
    })
    
  }) |> 
  ungroup()


greenCom_peak_areas <- GreenCom |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(GreenCom_peaks), function(i){
      
      peak <- GreenCom_peaks[i,]
      
      area <- riemann_integrate(
        .x,
        rt1 = peak$start_time,
        rt2 = peak$end_time
      )
      
      tibble(
        peak_number = peak$peak_group,
        start_time = peak$start_time,
        end_time = peak$end_time,
        apex_time = peak$apex_time,
        area = area
      )
      
    })
    
  }) |> 
  ungroup()

green_peak_areas <- bind_rows(greenCYP_peak_areas, greenUGT_peak_areas, greenCom_peak_areas)