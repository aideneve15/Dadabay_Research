"
Peak dataframes have the following information;
peak number, start time, end time as well as an apex_time

This script should take the start and end time for all peaks and 
integrate from the corresponding species.

Group by peak sample id and peak number and integrate
"

riemann_integrate <- function(df, rt1, rt2){
  # Pull numeric vectors directly
  time_vec <- as.numeric(df$time)
  intensity_vec <- as.numeric(df$intensity)
  
  # Filter indices
  keep <- time_vec >= rt1 & time_vec <= rt2
  time_vec <- time_vec[keep]
  intensity_vec <- intensity_vec[keep]
  
  # Need at least 2 points
  if(length(time_vec) < 2) return(0)
  
  dt <- diff(time_vec)
  intensity_mid <- (intensity_vec[-1] + intensity_vec[-length(intensity_vec)]) / 2
  sum(dt * intensity_mid, na.rm = TRUE)
}

allLong;

EarlyCYP <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP', species == 'early')

EarlyUGT <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'UGT', species == 'early')

EarlyCom <- allLong |> 
  mutate(sample_id = paste(sample, species, enzyme, time_point, sep = '_')) |> 
  dplyr::filter(enzyme == 'CYP+UGT', species == 'early')

EarlyCom_peaks;
EarlyCYP_peaks;
EarlyUGT_peaks;

library(purrr)
earlyCYP_peak_areas <- EarlyCYP |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(EarlyCYP_peaks), function(i){
      
      peak <- EarlyCYP_peaks[i,]
      
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

earlyUGT_peak_areas <- EarlyUGT |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(EarlyUGT_peaks), function(i){
      
      peak <- EarlyUGT_peaks[i,]
      
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


earlyCom_peak_areas <- EarlyCom |> 
  group_by(species, enzyme, time_point, sample) |> 
  group_modify(~{
    
    map_dfr(1:nrow(EarlyCom_peaks), function(i){
      
      peak <- EarlyCom_peaks[i,]
      
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

early_peak_areas <- bind_rows(earlyCYP_peak_areas, earlyUGT_peak_areas, earlyCom_peak_areas)