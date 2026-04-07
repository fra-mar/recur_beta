# Dataset creator forlater brms analysis
library('dplyr')

roc2microM = 1.887505
sug2microM = .499451
if (exists('adf')){rm(adf)}

dataSourceDir = '~/Code/recur/'    # Where folder simulatios is
pwd = '~/Code/recur/statAnalysis/' # Where to read/write mydata.csv

setwd(pwd)

# Check if data exists, otherwise create it---------------------------------
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
    print(summary(tdf_params))
    joined = left_join(tdf_params,tdf_tofAn,by='id')
    joined = left_join(joined, tdf_forReg, by='id')
    
    if (!(exists('mydata'))){
      print('starting mydata')
      mydata = joined
    } else {
      print('continuing mydata')
      mydata = tibble(rbind(mydata, joined))
    }

  } #---------- End of for loop
  # Convert cumROC, A1roc and A2roc to mg/Kg for
  # 1. correction by weight
  # 2. readers better intuition mg than microM
  mydata = mydata |>
    mutate(cumROC = cumROC /(BW*roc2microM),
           A1ro   = A1ro / (BW*roc2microM),
           A2ro   = A2ro / (BW*roc2microM))
  print(summary(mydata))
  setwd(pwd)
  write.csv(mydata, 'mydata.csv', row.names = FALSE)
  #------------ End of if, when no gathered data exists

  } else { 
  mydata = tibble(read.csv('mydata.csv'))
}  

# Set factors
mydata = mydata |>
  mutate(catSUG = factor(catSUG, levels=c('adjusted','fixed')),
         RAC = factor(RAC, levels=c('Asian','nonAsian')),  
         SEV = factor(SEV, levels=c('False','True')))

#--------------Support functions----------------------------------------------
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

