########################### Description #######################################

  ## This script compares counts of sedimentary structures across the 
  ## Phanerozoic using Fisher's exact test of independence. For groups of 
  ## n by 2 categorical variables, post hoc testing consists of pairwise Fisher 
  ## tests with a correction for multiple comparisons (see below). Pairwise 
  ## comparisons are rendered as a compact letter display.

  ## For this test, data from the Phanerozoic are pooled into four categories: 
  ## Cambrian-Devonian, Carboniferous-early Permian, mid- Permian to 
  ## mid-Jurassic, and Late-Mesozoic to Modern. Occurrences of sedimentary 
  ## structures (T for tepees and P for pisoids) are "successes" while counts 
  ## where these features are absent (NT and NP) are "failures"

###################### Load data and specify databases ########################

  ## Clear any variables from the global environment.

  remove(list= ls())

  ## Load the rcompanion library. Must be installed to run the code.
  ## 
  library(rcompanion)

########################### Tepee analysis ####################################
  
  ## Specifies the input data to be tested. These data are shown in Table 2 in 
  ## the text.

  Input =("
                  Time                   T        NT
          Cambrian_Devonian             51        29
          Carboniferous_Early_Permian    4        19
          mid_Permian_mid_Jurassic      39        15
          mid_Jurassic_pressent         15        70
          ")


  ## Convert the input table to matrix form to pass to the fisher.test function

  Matriz = as.matrix(read.table(textConnection(Input),
                                header=TRUE, 
                                row.names=1))


  ## Run Fisher's exact test for independence specifying the two-sided p value. 
  ## Fisher's exact test was chosen rather than the Chi-square test or G test 
  ## because the total counts are <1000; see the discussion on sample size in 
  ## McDonald (2014). Note that the line 'alternative = "two-sided"' is not 
  ## strictly necessary for a table larger than 2x2. 

  fisher.test(Matriz,
              alternative="two.sided")

  ## NOTE: The post-hoc testing in the next section should only be used if the 
  ## p-value given above is significant at the .05 confidence level. 
  
  ## Conduct pairwise Fisher tests among different time designations. 
  ## "Method" specifies the method for adjusting the p-values for multiple 
  ## comparisons. The Holm method is suggested as one alternative to
  ## the Bonferroni correction for controlling the family-wise error rate 
  ## (FWER). See the documentation for the p.adjust function in the R 
  ## documentation. 
   
  PT = pairwiseNominalIndependence(Matriz,
                                 
                                   fisher = TRUE,
                                   gtest  = FALSE,
                                   method = "holm",
                                   chisq  = FALSE)

  ## Display relationships among pairwise comparisons as a compact letter 
  ## display. Specifics on the algorithm and reference can be found in the R 
  ## documentation for the cld function.

  cldList(comparison = PT$Comparison,
          p.value    = PT$p.adj.Fisher,
          threshold  = 0.05)     

########################### Pisoid analysis ###################################
  
  ## Specifies the input data to be tested. These data are shown in Table 2 in 
  ## the text.
  
  Input =("
                  Time                  P        NP
          Cambrian_Devonian             6        74
          Carboniferous_Early_Permian   3        20
          mid_Permian_Early_Jurassic    19       35
          mid_Jurassic_pressent         4        81
          ")
  ## Convert the input table to matrix form to pass to the fisher.test function
  
  Matriz = as.matrix(read.table(textConnection(Input),
                                header=TRUE, 
                                row.names=1))

  ## Run Fisher's exact test for independence specifying the two-sided p value. 
  ## Fisher's exact test was chosen rather than the Chi-square test or G test 
  ## because the total counts are <1000; see the discussion on sample size in 
  ## McDonald (2014). Note that the line 'alternative = "two-sided"' is not 
  ## strictly necessary for a table larger than 2x2.

  fisher.test(Matriz,
              alternative="two.sided")
  
  ## NOTE: The post-hoc testing in the next section should only be used if the 
  ## p-value given above is significant at the .05 confidence level. 
  
  ## Conducts pairwise Fisher tests among different time designations. 
  ## "Method" specifies the method for adjusting the p-values for multiple 
  ## comparisons. The Holm method is suggested as one alternative to
  ## the Bonferroni correction for controlling the family-wise error rate 
  ## (FWER). See the documentation for the p.adjust function in the R 
  ## documentation.
  
  PT = pairwiseNominalIndependence(Matriz,
                                 
                                   fisher = TRUE,
                                   gtest  = FALSE,
                                   method = "holm",
                                   chisq  = FALSE)

  ## Display relationships among pairwise comparisons as a compact letter 
  ## display. Specifics on the algorithm and reference can be found in the R 
  ## documentation for the cld function.

  cldList(comparison = PT$Comparison,
          p.value    = PT$p.adj.Fisher,
          threshold  = 0.05)     
