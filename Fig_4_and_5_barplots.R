################# Barplots for tepee and pisoid data ########################### 

  ## This script plots total counts (Fig. 4a,b) and normalized counts 
  ## (Fig. 5a,b) for both tepees and pisoids

###################### Load data and specify databases #########################

  ## Loads input file with data for all time slices. Default is based on 
  ## Ridgwell (2005).

  ## Clear any variables from the global environment.

  remove(list= ls())

  ## Optional: Set the working directory. The working directory must include 
  ## the input file(s) (with headers) saved as a .csv file. 

  #setwd("............")
  
  ## Load libraries. All of these must be installed to run the code. 

  library(plyr)
  library(ggplot2)
  library(scales)  # for percentage scales

  ## Load data from the CSV file

  data <- read.csv("database_entries.csv")

  get_bins <-aggregate(data, list(bin = data$bin), FUN = mean)

  summary <- ddply(data, .(bin), summarize, 
                 P_yes = sum(Pisoids == 'Yes'),
                 P_no  = sum(Pisoids == 'No'),
                 T_yes = sum(Tepees  == 'Yes'),
                 T_no  = sum(Tepees  == 'No'))


  T_yes <- data.frame(summary$bin, summary$T_yes)
  summary$all <- summary$P_yes + summary$P_no
  summary$P_norm <- 100*summary$P_yes/summary$all
  summary$T_norm <- 100*summary$T_yes/summary$all

  to_plot <- data.frame(summary, get_bins$Lower, get_bins$Upper,
                      get_bins$Center, get_bins$Range)

################ Collect and plot data #####################################

  ## Raw counts for tepees (Fig. 4a)

  ggplot(to_plot, aes(ymin = 0)) +
  
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, ymax = T_yes),
                fill = 'blue', alpha = .6) + 
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, 
                ymin = T_yes, ymax = all),
                fill = 'grey') +
    
  scale_x_reverse() +
  theme_classic() 

  ## Raw counts for pisoids (Fig. 4b)

  ggplot(to_plot, aes(ymin = 0)) +
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, ymax = P_yes),
            fill = 'green', alpha = .6) + 
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, 
                ymin = P_yes, ymax = all),
            fill = 'grey') +
  scale_x_reverse() +
  theme_classic() 

  ## Normalized tepees (Fig. 5a)

  ggplot(to_plot, aes(ymin = 0)) +
    
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, ymax = T_norm),
          fill = 'blue', alpha = .6) + 
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, 
                ymin = T_norm, ymax = 100),
            fill = 'grey') +
  scale_x_reverse() +
  theme_classic() 

  ## Normalized pisoids (Fig. 5b)

  ggplot(to_plot, aes(ymin = 0)) +
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, ymax = P_norm),
            fill = 'green', alpha = .6) + 
  geom_rect(aes(xmin = get_bins.Upper, xmax = get_bins.Lower, 
                ymin = P_norm, ymax = 100),
            fill = 'grey') +
  scale_x_reverse() +
  theme_classic() 
