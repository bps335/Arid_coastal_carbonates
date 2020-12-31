######### Carbonate saturation states in modern coastal brines #################

  ## This script calculates aragonite saturation states in modern coastal brines
  ## surveyed from the literature and generates the violin plot in Fig. 2a.

###################### Load data and specify databases #########################

  ## Clear any variables from the global environment.

  remove(list= ls())

  ## Optional: Set the working directory. The working directory must include 
  ## the input file(s) (with headers) saved as a .csv file. 

  #setwd("............")

  ## Load libraries. All of these must be installed to run the code. 

  library(phreeqc)
  library(ggplot2)
  library(easyphreeqc)

  ## Specify the PHREEQC database file for the model and allow R to 
  ## query the output of the calculations. The default is the pitzer database 
  ## (pitzer.dat). 

  phrLoadDatabaseString(pitzer.dat)

  ## Load data from the CSV file

  data <- read.csv("modern_brines.csv")
  
  ## Save the total number of rows 
  
  length <-  dim(data)[1]  
  
  data$arag_omega <- rep(0,length)
  data$gyp_omega <- rep(0,length)
  
########################## Set PHREEQC parameters ##############################
  
  for (row in 1:length)
  {
    value = row
    skip_to_next <- FALSE
    
    input <- c(
            'SOLUTION 1              '                                         ,       
      paste('units                   ', data$units[value], sep = " ")          ,
      paste('pH                      ', data$pH[value], sep = " ")             ,
      paste('Alkalinity              ', data$Alkalinity[value], sep = " ")     ,
      paste('Ca                      ', data$Ca[value], sep = " ")             ,
      paste('Na                      ', data$Na[value], sep = " ")             ,
      paste('S(6)                    ', data$SO4[value], sep = " ")            ,
      paste('Cl                      ', data$Cl[value], sep = " ")             ,
      paste('Mg                      ', data$Mg[value], sep = " ")             ,
      paste('temp                    ', data$Temp[value], sep = " ")           ,  
      
      # Amend the PHREEQC database  
      
      '  SOLUTION_SPECIES       ', # Taken from the phreeqc.dat file 
      
      '  Ca+2 + CO3-2 = CaCO3   '                                              ,
      '  log_k	   3.224        '                                              ,
      '  delta_h   3.545   kcal '                                              ,
      
      '  SELECTED_OUTPUT'                                                      ,
      '  -si calcite'                                                          ,
      '  -si aragonite'                                                        ,
      '  -si gypsum'                                                           ,
      '  -si halite'                                                           ,
      '  -si natron'                                                           ,
      '  -ph  true'                                                            , 
      '  -alkalinity  true'                                                    ,
      '  -totals Ca Mg Na K S(6) Cl C(4)'                                      ,
      '  -water'                                                               )
    
    ## Run the first PHREEQC code and query the output.
    
    error <- tryCatch(phrRunString(input), error = function(e) e) 
    
    if ( length(error) == 0){
      
      phrRunString(input)  
      output <- phrGetSelectedOutput()
      

      data$arag_omega[value]   <-   10^(tail(output$n1$si_aragonite, n = 1))
      data$gyp_omega[value]   <-   10^(tail(output$n1$si_gypsum, n = 1))
    }
  }

########################## Cull results and plot ##############################

  ## Remove zeros (these are cells for which the program was unable to 
  ## converge on an answer)

  data<-subset(data, data$arag_omega>0) 
  
  ## Focus on brines below gypsum saturation
  
  data<-subset(data, data$gyp_omega<1) 
  

  #data<-subset(data, data$omega_arag<20)
  
  
  ggplot (data, aes(x = Sample.depth, 
                    y = log10(arag_omega), 
                    fill = Sample.depth)) +
    
    annotate(geom = "rect", 
             ymin = log10(3), xmin = -Inf,
             ymax = log10(5), xmax = Inf,
             fill = "dark blue", alpha = 0.4)  +
    
    geom_violin(trim = FALSE) + geom_boxplot(width=0.1, fill = "white") 
  
  
  
############################# Summary statistics ###############################
  
  surface   <-subset(data, data$Sample.depth == "surface")
  porewater <-subset(data, data$Sample.depth == "porewater")
  
  surf_median = median(surface$arag_omega)
  pore_median = median(porewater$arag_omega)
  
  ## Perform t-test w/unequal variances on calculated saturation states. 
  
  t.test(log10(porewater$arag_omega),log10(surface$arag_omega))