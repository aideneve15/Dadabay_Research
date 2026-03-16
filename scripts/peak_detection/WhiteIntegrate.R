"
Peak dataframes have the following information;
peak number, start time, end time as well as an apex_time

This script should take the start and end time for all peaks and 
integrate from the corresponding species.

Group by peak sample id and peak number and integrate
"

allLong;

WhiteCYP <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP', species == 'white')

WhiteUGT <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'UGT', species == 'white')

WhiteCom <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP+UGT', species == 'white')

WhiteCom_peaks;
WhiteCYP_peaks;
WhiteUGT_peaks;

library(purrr)
whiteCYP_peak_areas <- WhiteCYP |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WhiteCYP_peaks), function(i){
      
      peak <- WhiteCYP_peaks[i,]
      
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

whiteUGT_peak_areas <- WhiteUGT |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WhiteUGT_peaks), function(i){
      
      peak <- WhiteUGT_peaks[i,]
      
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


whiteCom_peak_areas <- WhiteCom |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(WhiteCom_peaks), function(i){
      
      peak <- WhiteCom_peaks[i,]
      
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

white_peak_areas <- bind_rows(whiteCYP_peak_areas, whiteUGT_peak_areas, whiteCom_peak_areas)