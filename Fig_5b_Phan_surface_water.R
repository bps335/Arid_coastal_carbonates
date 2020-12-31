########### Phanerozoic trend for surface water-kinetic  model ################# 

  ## This script reproduces the model results in Fig. 5b.

###################### Load data and specify databases #########################

  ## Loads input file with data for all time slices. Default is based on 
  ## Ridgwell (2005).

  ## Clear any variables from the global environment.

  remove(list= ls())

  ## Optional: Set the working directory. The working directory must include 
  ## the input file(s) (with headers) saved as a .csv file. 

  setwd("F:/New Backup/PhD Work/Tepee_GRL/Github")
  getwd()

  ## Load libraries. Both of these must be installed to run the code. 

  library(phreeqc)
  library(ggplot2)

  ## Specify the PHREEQC database file for the model and allow R to 
  ## query the output of the calculations. The default is the pitzer database 
  ## (pitzer.dat). 

  phrLoadDatabaseString(pitzer.dat)

  ## Load data from the CSV file

  data <- read.csv("model_input.csv")

  ## Unit conversion for DIC and alkalinity data (reported in micromoles)
  ## to millimoles 

  data$DIC_surf <- data$DIC_surf/1000
  data$Carb_alk <- data$Carb_alk/1000

###################### Set model parameters ##############################

  ## Set the total number of volumes flushed through the system. This 
  ## paramters non-dimensionalized by normalizing the total volumes, V, to 
  ## the starting volume vO. This approach is similar to Sanford and Wood 
  ## (1990) and Jacobsen and Kaufman (1999); using dimensionless volume 
  ## abstracts the time-dependant aspects of the system.

  vmax = 100                

  ## Set the step increment. To ensure convergance to the steady state, 
  ## compare the model output to the analytic solution for conservative 
  ## elements (e.g., Cl) from Sanford and Wood (1990).

  dv = 2     
  
  ## Set the number of steps in the model based on the maximum number of 
  ## volumes and the step size.
    
  nv = vmax/dv            
  
  ## Set the paramter f, which partitions outgoing water fluxes between 
  ## evaporation and outflow. Higher f values specify more outflow (f*dv) 
  ## relative to evaporation ((1-f)*dv). 

  f = .75 
  
  ## Note that if we assume that carbonate precipitation has a negligable 
  ## effect on salinity, we can estimate the steady-state salinity using 
  ## the analytic solution for conservative solutes given by Sanford and 
  ## Wood (1990). This allows us to relate the partitioning of outgoing 
  ## fluxes to the salinity in the steady-state solution.
  
  salinity = 35/f
  
############################# Run the model ##############################

  ## Set the temperature loop and specify the lower, middle, and upper 
  ## temperature bounds(Celcius) for the model in degrees Celcius
  
  t_low = 30
  t_mid = 40
  t_high = 50
  
  for (t in c(t_low, t_mid, t_high))
  {

  ## Query the number of rows in the input file and pre-allocate a vector 
  ## for the steady-state calcite saturation state. 
    
  n = dim(data)[1]  
  cal_out = rep(0,n)
  
  ## Set the age loop.
  
  for (time in 1:n)
  {
  bin = time
  
  ## Set the initial seawater composition for the designated row.
  ## Since the carbonate system has two degrees of freedom, two of the
  ## six parameters (DIC, alkalinity, pH, CO2, CO3, and HCO3) must be
  ## specified (Zeebe and Wolf-Gladrow, 2001). The two variables used here
  ## are alkalinity and carbonate ion (CO3) concentration. The convention in 
  ## PHREEQC is to specify a 'dummy' pH (e.g., pH = 8) and then specify the 
  ## saturation state of a mineral (calcite) on the same line. Since Ca is also
  ## specified, this sets the carbonate ion concentration. 
  
  ## I prefer this approach because it avoids conversion issues with pH scales.
  ## PHREEQC apparently defines pH by the activity of hydrogen, which will 
  ## yield inconsistent results if the input pH is on the seawater scale or the
  ## total scale. It is then convenient to use DIC and alkalinity alkalinity 
  ## because they are conservative quantities that work well with the MIX 
  ## block. The starting temperature, 15.7 C, is the globally-averaged value 
  ## in the model from Ridgwell (2005).    
  
  input <- c(
    '  SOLUTION 1            '                                              ,       
    '  units   mmol/kgw      '                                              ,
    paste ('C(4)             ', data$DIC_surf[bin], sep = ' ')              ,
    paste ('pH    8  calcite ', data$SI_cal[bin],sep = ' ')                 ,
    paste ('Ca               ', data$Ca[bin], sep = ' ')                    ,
    paste ('Mg               ', data$Mg[bin], sep = ' ')                    ,
    paste ('Na               ', data$Na[bin], sep = ' ')                    ,
    paste ('S(6)             ', data$SO4[bin], sep = ' ')                   , 
    paste ('K                ', data$K[bin], sep = ' ')                     ,
    paste ('Cl               ', data$Cl[bin], sep = ' ')                    ,
           'temp         15.7'                                              ,
    
    # Amend the PHREEQC database  
    
    ' SOLUTION_SPECIES       ', # Taken from the phreeqc.dat file 
    
    ' Ca+2 + CO3-2 = CaCO3   '                                              ,
    ' log_k	   3.224         '                                              ,
    ' delta_h   3.545   kcal '                                              ,
   
  ## Set the parameter to send to the output file. There are several 'extra'
  ## outputs that are not critical to this project, but are implemented as a
  ## seed for future projects. 
  
    '  SELECTED_OUTPUT'                                                     ,
    '  -high_precision true'                                                ,
    '  -si calcite'                                                         ,
    '  -si aragonite'                                                       ,
    '  -si gypsum'                                                          ,
    '  -si halite'                                                          ,
    '  -si natron'                                                          ,
    '  -ph  true'                                                           , 
    '  -alkalinity  true'                                                   ,
    '  -totals Ca Mg Na K S(6) Cl C(4)'                                     ,
    '  -water'                                                              )
  
  ## Run the first PHREEQC code and query the output.
  
  phrRunString(input)
  output <- phrGetSelectedOutput()
  
  ## Overwrite the carbonate alkalinity from the input file and replace it 
  ## with the alkalinity adjusted to achieve the desired saturation state.
  
  data$Carb_alk[bin] <- output$n1$Alk.eq.kgw*1000
  
  ## Copy the data into a dummy file. This file will be overwritten to 
  ## store the composition of evaporated seawater during model run. 
  ## Copying the data simply gives the first timestep the composition 
  ## of unaltered seawater. 
  
  data_evap <- data
  
  ## Set the loop for the open-system model. 
  for (i in 1:nv)  
  {   
  temp = t
  
  ## Set the amount of fluid in the control volume after the last evaporation
  ## step (Vf1) and the inflow of fresh seawater to replace volume lost during
  ## evaporation (Vf2).
  
  Vf1 = 1-f*dv
  Vf2 = 1*dv
  
  ## Sets the amount of water (in moles) to be removed from the system 
  ## during the evaporation phase. The amount is based on a control volume
  ## of 1 kg water, which is the default in PHREEQC. 
  
  C1 = (1-f)*dv*55.506
  
  ## Set the PHREEQC input. The first solution is the evaporated seawater 
  ## already present in the control volume. The second solution is the influx
  ## of regular, unevaporated seawater into the system. 

  input <- c(
    '  SOLUTION 1             '                                             ,
    '  units   mmol/kgw       '                                             ,
    paste ('C(4)              ', data_evap$DIC_surf[bin], sep = ' ')        ,
    paste ('Alkalinity        ', data_evap$Carb_alk[bin],sep = ' ')         ,
    paste ('Ca                ', data_evap$Ca[bin], sep = ' ')              ,
    paste ('Mg                ', data_evap$Mg[bin], sep = ' ')              ,
    paste ('Na                ', data_evap$Na[bin], sep = ' ')              ,
    paste ('S(6)              ', data_evap$SO4[bin], sep = ' ')             , 
    paste ('K                 ', data_evap$K[bin], sep = ' ')               ,
    paste ('Cl                ', data_evap$Cl[bin], sep = ' ')              ,
    paste ('temp              ', temp, sep = ' ')                           ,
   
     
    '  SOLUTION 2             '                                             ,
    '  units   mmol/kgw       '                                             ,
    paste ('C(4)              ', data$DIC_surf[bin], sep = ' ')             ,
    paste ('Alkalinity        ', data$Carb_alk[bin],sep = ' ')              ,
    paste ('Ca                ', data$Ca[bin], sep = ' ')                   ,
    paste ('Mg                ', data$Mg[bin], sep = ' ')                   ,
    paste ('Na                ', data$Na[bin], sep = ' ')                   ,
    paste ('S(6)              ', data$SO4[bin], sep = ' ')                  ,
    paste ('K                 ', data$K[bin], sep = ' ')                    ,
    paste ('Cl                ', data$Cl[bin], sep = ' ')                   ,
    paste ('temp              ', temp, sep = ' ')                           ,
    
  ## Mix the evaporated water (Vf1) with the appropriate volume of fresh 
  ## seawater (Vf2)
  
    '  Mix 1                  '                                             ,
    paste('1                  ', Vf1, sep = ' ')                            ,
    paste('2                  ', Vf2, sep = ' ')                            ,
    
  ## Remove water during the evaporation step
  
    '  REACTION                '                                            ,
    '  H2O     -1.0            '                                            ,
    paste(C1, ' moles in 1 step', sep = ' ')                                ,
  
  # Amend the PHREEQC database  
  
  ' SOLUTION_SPECIES       ', # Taken from the phreeqc.dat file 
  
  ' Ca+2 + CO3-2 = CaCO3   '                                              ,
  ' log_k	   3.224         '                                              ,
  ' delta_h   3.545   kcal '                                              ,
    

  ## Set the parameter to send to the output file.
   
    '  SELECTED_OUTPUT'                                                     ,
    '  -high_precision true'                                                ,
    '  -si calcite'                                                         ,
    '  -si aragonite'                                                       ,
    '  -si gypsum'                                                          ,
    '  -si halite'                                                          ,
    '  -si natron'                                                          ,
    '  -equilibrium_phases calcite gypsum halite'                           ,
    '  -ph  true'                                                           , 
    '  -alkalinity  true'                                                   ,
    '  -totals Ca Mg Na K S(6) Cl C(4)'                                     ,
    '  -water'                                                              )
  
  ## Run the PHREEQC step and query the output.
  
  phrRunString(input)
  output <- phrGetSelectedOutput()
  
  ## Update the composition of the control volume to feed to the next step in
  ## the loop.

  data_evap$Ca[bin] <- tail(output$n1$Ca.mol.kgw.*1000, n = 1)
  data_evap$Na[bin] <- tail(output$n1$Na.mol.kgw.*1000, n = 1)
  data_evap$Cl[bin] <- tail(output$n1$Cl.mol.kgw.*1000, n = 1)
  data_evap$Mg[bin] <- tail(output$n1$Mg.mol.kgw.*1000, n = 1)
  data_evap$SO4[bin] <- tail(output$n1$S.6..mol.kgw.*1000, n =1)
  data_evap$Carb_alk[bin] <- tail(output$n1$Alk.eq.kgw.*1000, n =1)
  data_evap$DIC_surf[bin] <- tail(output$n1$C.4..mol.kgw.*1000, n =1)
  data_evap$K[bin] <- tail(output$n1$K.mol.kgw.*1000, n =1)
  
  ## Close the loop for the steady-state model.
  }
  
  ## Query the final (steady-state) saturation state from the loop and store
  ## it in a pre-allocated vector. 
  
  cal_out[bin] <- tail(output$n1$si_calcite, n = 1)
  
  ## Close the age loop for a single temperature. 
  }
  
  ## Convert saturation indices to saturation states and store them in a 
  ## vector named by temperature
  
  if (t == t_low){
    assign('phan_low',(10^cal_out))
  } 
  
  else if (t == t_mid){
    assign('phan_mid',(10^cal_out))  
  }
  else if (t == t_high){
    assign('phan_high',(10^cal_out))  
  }
  
  ## Close the temperature loop.
  }  
  
################ Collect and plot data #####################################

  ## Collect saturation state outputs, convert them to 'normalized 
  ## precipitation rates' (see SI), then add them to the data frame.
   
  phan_low_out = (phan_low-1)^3.12/(phan_mid[n]-1)^3.12
  phan_mid_out = (phan_mid-1)^3.12/(phan_mid[n]-1)^3.12
  phan_high_out = (phan_high-1)^3.12/(phan_mid[n]-1)^3.12

  data<-data.frame(data,phan_low_out,phan_mid_out,phan_high_out)
  
  ## Plot the outputs, starting with the average temperature case.
  
  ggplot(data, aes(x = Age.bin, y = phan_mid_out)) + 
    geom_line(size = 1.2, color = 'black') +
    
  ## Plot the upper and lower temperature bounds and fill the area between
  ## as an estimate of the uncertainty. 
  
    geom_ribbon(aes(ymin=phan_low_out, ymax = phan_high_out), 
                fill = 'blue', alpha=0.3) +
    scale_x_reverse() +
    theme_classic() +
  
  ## Add labels and tidy up the plot for export.
  
    labs (x = 'Age (Ma)', y = 'Precipitation factor', 
          title = 'Surface water model') +
    theme(plot.title = element_text(size = rel(1.5), hjust = .5)) + 
    theme(axis.title  = element_text(size = rel(1.2), hjust = .5)) 
  
