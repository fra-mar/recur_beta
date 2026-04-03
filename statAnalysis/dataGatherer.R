# Dataset creator forlater brms analysis
library('dplyr')

roc2microM = 1.887505
sug2microM = .499451
if (exists('adf')){rm(adf)}

# Check if data exists, otherwise create it
dataSourceDir = '~/Code/recur/'


if (!('mydata.csv' %in% list.files())) {
  
  print('Datafile doesnt exist. Creating it')
  
  setwd(dataSourceDir)
  #scenarios = tibble(read.csv('scenarios.csv'))
  
  for (scenario in list.dirs('simulations', recursive=FALSE)){
   print(paste0('Working on ',scenario))
    # Read tofAnalyzer
    tdf_tofAn  = tibble(
      read.csv(
        paste0(scenario,'/analyzedTOF.csv'))) |>
      mutate(id = as.integer(substr(sim.id,3,8))) |>
      select(id, spntRecover,isRNB, t90, tIniRNB, durRNB, lowestTOF )  
 
    # Read forRegression
    tdf_forReg = tibble(
      read.csv(
        paste0(scenario,'/forRegression.csv'))) |>
      mutate(id = as.integer(id)) |>
      rename(c(cumROC='dROC'))
  
    # Read params
    tdf_params  = tibble(
      read.csv(
        paste0(scenario,'/params.csv'))) |>
      select(!simId)
    
    joined = left_join(tdf_params,tdf_tofAn,by='id')
    joined = left_join(joined, tdf_forReg, by='id') |>
      mutate(cumROC = cumROC/BW) #cummulative ROC(microM) per Kg BW
    
    if (!(exists('adf'))){
      print('starting adf')
      adf = joined
    } else {
      print('continuing adf')
      mydata = tibble(rbind(adf, joined))
    }
    
  
  } #---------- End of for loop
  
  mydata = mydata |>
    mutate(rSUG_remROC = dSUG/(A1ro+A2ro),  # given SUG/remainingROC
           A1ro = A1ro / roc2microM, #ROC1 -> mg, so reader understands)
           A2ro = A2ro / roc2microM)
  
  write.csv(mydata, 'mydata.csv', row.names = FALSE)
  #------------ End of if, when no gathered data exists

  } else { 
  mydata = tibble(read.csv('mydata.csv'))
}  

# Set factors
mydata = mydata |>
  mutate(catSUG = factor(catSUG, levels=c('adjusted','fixed')),
         RAC = factor(RAC),
         SEV = factor(SEV))

#--------------Support functions------------------
qtMe = function(x, probs = c(0.05, 0.5, 0.95)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('p5','p50','p95')
  )
}

mySummary = function(df) {
  summarized = df |> summarise(across(where(is.numeric), .fns = 
                          list(median = ~quantile(.,.5),
                               p5 = ~quantile(., 0.05),
                               p95 = ~quantile(., 0.95)))) |>
  tidyr::pivot_longer(everything(), 
                      names_sep='_', 
                      names_to=c('variable', '.value'))
  
  return (summarized)
}
mySummary2 = function(df) {
  summarized = df |> summarise(across(where(is.numeric), .fns = 
                                        list(median = ~quantile(.,.5),
                                             p5 = ~quantile(., 0.05),
                                             p95 = ~quantile(., 0.95)))) |>
    tidyr::pivot_longer(everything())
  
  return (summarized)
}

