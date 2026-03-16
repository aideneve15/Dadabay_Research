"
Peak dataframes have the following information;
peak number, start time, end time as well as an apex_time

This script should take the start and end time for all peaks and 
integrate from the corresponding species.

Group by peak sample id and peak number and integrate
"

allLong;

WyomingCYP <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP', species == 'wyoming')

WyomingUGT <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'UGT', species == 'wyoming')

WyomingCom <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP+UGT', species == 'wyoming')

WyomingCom_peaks;
WyomingCYP_peaks;
WyomingUGT_peaks;

library(purrr)
wyomingCYP_peak_areas <- WyomingCYP |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WyomingCYP_peaks), function(i){
      
      peak <- WyomingCYP_peaks[i,]
      
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

wyomingUGT_peak_areas <- WyomingUGT |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WyomingUGT_peaks), function(i){
      
      peak <- WyomingUGT_peaks[i,]
      
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


wyomingCom_peak_areas <- WyomingCom |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WyomingCom_peaks), function(i){
      
      peak <- WyomingCom_peaks[i,]
      
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

wyoming_peak_areas <- bind_rows(wyomingCYP_peak_areas, wyomingUGT_peak_areas, wyomingCom_peak_areas)