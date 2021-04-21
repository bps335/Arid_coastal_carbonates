###################### DIC-alkalinity contour plot ############################# 

  ## This script creates a contour polot of calcite Omega as a function of DIC 
  ## and alkalinity for Ca = 10 mmol as shown in Fig. 3c.  

###################### Load data and specify databases #########################

  ## Loads input file with data for all time slices. Default is based on 
  ## Ridgwell (2005).

  ## Clear any variables from the global environment.

  remove(list= ls())

  ## Optional: Set the working directory. The working directory must include 
  ## the input file(s) (with headers) saved as a .csv file. 

  #setwd("............")
  setwd("F:/New Backup/PhD Work/Tepee_GRL/Revisions")
  
  ## Load libraries. All of these must be installed to run the code.
  
  library(easyphreeqc)
  library(plotly)
  library(colorRamps)
  library(processx)
  library(orca)

  ## Load data from the CSV file

  data <- read.csv("model_input.csv")

########################## Set plot parameters #################################

  ## create data.frame with arguments to solution()

  m = 75
  n = 75

  DIC_min = 1.0         ##[mmol/kgw]
  DIC_max = 8.5         ##[mmol/kgw]

  Alk_min = 1.0         ##[mmol/kgw]
  Alk_max = 8.5         ##[mmol/kgw]

  DIC <- seq(from = DIC_min, to = DIC_max, by = (DIC_max - DIC_min)/(m - 1))
  Alk <- seq(from = Alk_min, to = Alk_max, by = (Alk_max - Alk_min)/(n - 1))

  grid <- expand.grid(DIC = DIC, Alk = Alk)
  end = m*n

  ## pass plot parameters to solution_info
  
  solution_info <- data.frame(
  units = 'mmol/kgw',
  "C(4)" =  grid$DIC,
  "S(6)" = 28,
  Na = 469,
  Mg = 53,
  K = 11, 
  Ca = 10.7,
  Cl = 546,
  Alkalinity = grid$Alk,
  temp = 25,
  .number = seq(from = 1, to =end, by =1), check.names = FALSE)

  ## use mlply to turn the grid into a list()
  
  solution_list <- plyr::mlply(solution_info, solution)

  ## pass to phreeqc()
  output <-phreeqc(
  solution_list,
  selected_output(molalities = c("CO3-2", "CO2", "HCO3-"),
                  totals = "C(4)", si = "calcite",
                  pH = TRUE, alkalinity = TRUE),
  db = "pitzer"
  )
  
  ## Convert saturation indices to saturation states and store them in a 
  ## vector 
  omega_calcite<- 10^output$si_calcite
  
  ## Reshape vector into grid with the same size as the input data

  z <- t(matrix(omega_calcite, m, n))
  
#################### Collect and plot data #####################################

  ## Plot results using the plot_ly library

  x <- list(title = "DIC (mmol)")
  y <- list(title = "Total alkalinity (mmol)")
  title<- list(text = "Omega_calcite contours for DIC and TA", x = .55)

  colors <- rev(matlab.like(100))

  plot <- plot_ly()

  plot <- plot %>% add_trace(x = DIC, y = Alk,
                      z = z, type = "contour", line = list(smoothing = .85),
                contours = list(start = -1, end = 40, size = .2,
                                showlines = FALSE, opacity = 0.6))


  plot <- plot %>% add_trace(x = DIC, y = Alk,
                        z = z, type = "contour", contours = list( start = 1,
                        end = 1, size = 1,  showlabels = FALSE, 
                        size = 3, coloring = "none"),
                        line = list(smoothing = 0.85, color = "red", width = 2)) 

  plot <- plot %>% add_trace(x = DIC, y = Alk,
                           z = z, type = "contour", contours = list( start = 5,
                           end = 20, size = 5,  showlabels = FALSE, 
                           size = 3, coloring = "none"),
                           line = list(smoothing = 0.85, color = "white", 
                                                         width = 1)) 


  plot <- plot %>% add_trace( data = data, x = data$DIC_surf/1000,
                                 y = data$Carb_alk/1000, 
                                type = "scatter", mode = "markers", 
                              symbol = "x", marker = list(color = "white"))

  plot


  a = .2
  plot %>% hide_colorbar() %>% hide_legend() %>%
  layout(xaxis = x, yaxis = y, title = title )

