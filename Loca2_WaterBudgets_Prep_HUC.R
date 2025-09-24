
# Libraries

# Tidyverse resources

library(package = "tidyverse") # Multiple Tidyverse Resources
library(package = "lubridate") # Date-Time Control


library(package = "ClimClass")



# User Modification Area

# Climate Division Selection

# Spatial Aggregation Statistic
#   MEAN, P000, P025, P050, P075, P100

target_aggregation_statistic = "MEAN"


# Pull the Metadata for the Climate Divisions

huc_metadata_file = str_c("https://thredds.ias.sdsmt.edu:8443/thredds/fileServer/" , 
                          "LOCA2/" ,
                          "Specific_Regional_Aggregate_Sets/USGS_HUC08_Basins/" ,
                          "USGS_HUC08_LUT.RData", 
                          sep = "")

load(url(huc_metadata_file),
     verbose = TRUE)
closeAllConnections()

HUC_table_avail = read_csv("https://thredds.ias.sdsmt.edu:8443/thredds/fileServer/LOCA2/Specific_Regional_Aggregate_Sets/USGS_HUC08_Basins/HUC_table_avail.csv")

# Pull the huc data

HUC_table_avail = HUC_table_avail |> filter(Enhanced==1)

all_hucs = unique(HUC_table_avail$`HUC-08`)

remove(huc_metadata_file)

print(all_hucs)



for (target_huc in all_hucs) {
  
  loca_metadata = USGS_HUC08_LUT |>
    filter(huc08 == target_huc)
  
  
  
  target_climate_name = str_c(unique(loca_metadata$huc08), 
                              " HUC, ", 
                              unique(loca_metadata$huc08_name),
                              " (",
                              unique(loca_metadata$huc04_name),
                              ")", 
                              sep = "")
  
  print(target_climate_name)
  
  # Crack our LOCA File for our climate zone
  
  
  loca_file_name = str_c("https://thredds.ias.sdsmt.edu:8443/thredds/fileServer/" , 
                         "LOCA2/" ,
                         "Specific_Regional_Aggregate_Sets/USGS_HUC08_Basins/" ,
                         "R_Monthly_Files/LOCA2_V1_HUC08_MONTHLY_" , 
                         target_huc ,
                         ".RData", 
                         sep = "")
  
  load(url(loca_file_name),
       verbose = TRUE)
  
  loca2_monthly$Division = target_huc 
  
  loca2_monthly$Division = factor(x      = loca2_monthly$Division,
                                  levels = USGS_HUC08_LUT$huc08)
  closeAllConnections()
  
  
  
  
  loca2_monthly = loca2_monthly |> 
    filter(Percentile == target_aggregation_statistic) |>
    mutate("tasavg" = (tasmax + tasmin)/2)  |>
    rename("HUC08"  = "Division")     |>
    select("Scenario",
           "Model",
           "Member",
           "Percentile",
           "Year",
           "Month",
           "HUC08",
           "Time",
           "tasmax",
           "tasavg",
           "tasmin",
           "pr")
  
  # Get Model Values
  
  Models    = unique(loca2_monthly$Model)
  Scenarios = unique(loca2_monthly$Scenario)
  
  remove(loca_file_name)
  
  


  

  
  
  # Get Model Values
  
  Models    = unique(loca2_monthly$Model)
  Scenarios = unique(loca2_monthly$Scenario)
  
  
  
  # ClimClass Time Series
  first = TRUE
  
  for (model in Models) {
    
    local_loca2_monthly = loca2_monthly |>
      filter(Model == model)
    
    
    Local_Scenarios = unique(local_loca2_monthly$Scenario)
    print(str_c("   - ",model))
    
    for (scenario in Local_Scenarios[2:length(Local_Scenarios)]) { 
      print(str_c("     - ",scenario))
      thorntwaite_inputs = local_loca2_monthly |>
        filter((Scenario == "Historical") | (Scenario == scenario)) |>
        arrange(Time) |>
        rename("year"  = Year,
               "month" = Month,
               "Tn"    = tasmin,
               "Tm"    = tasavg,
               "Tx"    = tasmax,
               "P"     = pr)
      thorntwaite_budget_raw = thornthwaite(series   = thorntwaite_inputs,
                                            latitude = loca_metadata$huc08_mean_latitude,
                                            TAW      = loca_metadata$huc08_mean_mass_content_of_water_in_soil,
                                            snow.init = 0.)
      
      # Precipitation
      water_budget = t( as_tibble(thorntwaite_budget_raw$W_balance$Precipitation) )
      colnames(water_budget) = str_c(1:12)
      water_budget = water_budget %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      water_budget$Variable = "Precipitation"
      
      #Et0
      raw_wb2 = t( as_tibble(thorntwaite_budget_raw$W_balance$Et0) )
      colnames(raw_wb2) = str_c(1:12)
      raw_wb2 = raw_wb2 %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      raw_wb2$Variable = "Potential_Evap"
      
      water_budget = rbind(water_budget,raw_wb2)
      remove(raw_wb2)
      
      #Storage
      raw_wb2           = t( as_tibble(thorntwaite_budget_raw$W_balance$Storage) )
      colnames(raw_wb2) = str_c(1:12)
      raw_wb2 = raw_wb2 %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      raw_wb2$Variable = "Storage"
      
      water_budget = rbind(water_budget,raw_wb2)  
      remove(raw_wb2)
      
      #'Prec. - PotEvap.'
      raw_wb2 = t( as_tibble(thorntwaite_budget_raw$W_balance$'Prec. - Evap.') )
      colnames(raw_wb2) = str_c(1:12)
      raw_wb2 = raw_wb2 %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      raw_wb2$Variable = "Prec_m_PE"
      
      water_budget = rbind(water_budget,raw_wb2)  
      remove(raw_wb2)
      
      #Deficit
      raw_wb2 = t( as_tibble(thorntwaite_budget_raw$W_balance$Deficit) )
      colnames(raw_wb2) = str_c(1:12)
      raw_wb2 = raw_wb2 %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      raw_wb2$Variable = "Deficit"
      
      water_budget = rbind(water_budget,raw_wb2)
      remove(raw_wb2)
      
      #Surplus
      raw_wb2 = t( as_tibble(thorntwaite_budget_raw$W_balance$Surplus) )
      colnames(raw_wb2) = str_c(1:12)
      raw_wb2 = raw_wb2 %>%
        as.data.frame %>% 
        rownames_to_column(.,
                           var = 'year')
      raw_wb2$Variable = "Surplus"
      
      water_budget = rbind(water_budget, raw_wb2)
      remove(raw_wb2)
      
      water_budget = gather(data  = water_budget,
                            key   = month,
                            value = "value",
                            str_c(1:12))
      
      water_budget$Date = as.Date(str_c(water_budget$year,
                                        "-",
                                        water_budget$month,
                                        "-15",
                                        sep = ""))
      
      water_budget = spread(data = water_budget,
                            key  = "Variable",
                            value = "value")
      
      water_budget = water_budget %>% arrange(Date)
      
      # finish the budget by critical parameters
      
      # calculate evapotransporation
      water_budget = water_budget %>% 
        mutate(Evaporation = Potential_Evap - Deficit)
      
      # calculate precipitation - true evaporation
      water_budget = water_budget %>% 
        mutate(Prec_m_Evap = Precipitation - Evaporation)
      
      
      
      # calculate recharge by calculating the increase in soil storage from rainfall
      water_budget = water_budget %>% 
        mutate(Recharge = c(NA, diff(x = Storage, 
                                     lag = 1))) %>% 
        mutate(Recharge = ifelse(test = Recharge>0, 
                                 yes  = Recharge, 
                                 no   = 0))
      
      # separate recharge from surplus in teh water budget
      water_budget = water_budget %>% 
        mutate(Surplus = Surplus - Recharge)  %>% 
        mutate(Surplus = ifelse(test = Surplus>0, 
                                yes  = Surplus, 
                                no   = 0))
      
      # calculate recharge by calculating the increase in soil storage from rainfall
      water_budget = water_budget %>% 
        mutate(Snowpack = Prec_m_Evap - Recharge - Surplus)  %>% 
        mutate(Snowpack = ifelse(test = Snowpack>0, 
                                 yes  = Snowpack, 
                                 no   = 0))  
      
      # repair precip-PE (seems to be a typo in the original ClimClass Code)
      water_budget = water_budget %>% 
        mutate(Prec_m_PE = Precipitation - Potential_Evap)    
      
      
      water_budget$Temp_Avg = round(thorntwaite_inputs$Tm[1:length(water_budget$Precipitation)],2)
      water_budget$Temp_Max = round(thorntwaite_inputs$Tx[1:length(water_budget$Precipitation)],2)
      water_budget$Temp_Min = round(thorntwaite_inputs$Tn[1:length(water_budget$Precipitation)],2)
      
      water_budget$Time             = thorntwaite_inputs$Time
      water_budget$Snowpack         = round(water_budget$Snowpack,5)
      water_budget$Scenario         = thorntwaite_inputs$Scenario
      water_budget$Model            = thorntwaite_inputs$Model
      water_budget$Member           = thorntwaite_inputs$Member
      water_budget$Percentile       = thorntwaite_inputs$Percentile
      water_budget$HUC08            = thorntwaite_inputs$HUC08
      
      
      
      water_budget = water_budget %>% select(Scenario,
                                             Model,
                                             Member,
                                             Percentile,
                                             Time,
                                             HUC08,
                                             Temp_Min,
                                             Temp_Avg,
                                             Temp_Max,
                                             Precipitation,
                                             Potential_Evap,
                                             Evaporation,
                                             Deficit,
                                             Storage,
                                             Snowpack,
                                             Recharge,
                                             Surplus,
                                             Prec_m_PE,
                                             Prec_m_Evap)
      
      if (scenario != Local_Scenarios[2]) {
        water_budget = water_budget |>
          filter(Scenario != "Historical")
      }
      
      if (first) {
        LOCA2_Water_Budget = water_budget
        first = FALSE
      } else {
        LOCA2_Water_Budget = rbind(LOCA2_Water_Budget,
                                   water_budget)
        
      }
    }
    
  }
  
  
  
  loca_file_name = str_c("./LOCA2_V1_LOCA2_THORNTHWAITE_" , 
                         target_huc ,
                         ".RData", 
                         sep = "")
  
  
  save(LOCA2_Water_Budget, file = loca_file_name)
  closeAllConnections()
  
}
