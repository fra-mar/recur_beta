# Goal1: Summarize results regarding TOF and pkpd parameters, for each scenario
rm(list=ls())
library(dplyr)
library(ggplot2)
#library(patchwork)

# Supporting functions---------------------------------------------------------
quantileMe = function(sr, srName){
  toReturn = quantile(sr, c(0.05,.5,.95))
  names(toReturn) = paste0(srName,'_',c('p5','p50','p95'))
  return (toReturn)}
  
whichEvent = function(spntRecover, isRNB, isDIP, didRecover){
    if (spntRecover==TRUE){event = 1}
    else if (isRNB==TRUE){event=2}
    else if (isDIP==TRUE){event=3}
    else if (didRecover==FALSE){event=4}
    else {event=0} 
  return (event)
}


# Functions for human readable data
npMe = function(n, p){
  toRet = paste0(as.character(round(n,3)),
                 ' (',
                 as.character(round(p*100,1)),
                 '%)')
  return (toRet)
}

qMe = function(p5,p50,p95){
  if (any(is.na(c(p5,p50,p95)))) {
    toRet = '-'
  }
  else{
    toRet = paste0(as.character(round(p50,1)),
                   ' (',
                   as.character(round(p5,1)),
                   ',',
                   as.character(round(p95,1)),
                   ')')
  }
  return(toRet)
}
# For these params  variability was reported.
targetPars = c("simId", "CLro" ,"V1ro", "V2ro", "CLsu", 
               "ke0", "EC50", "Hill", "ks")
# To visuzalize in last table scenarios with no RNB events
scenWithNoEvents = c()

scenDirs = list.dirs('simulations', recursive = FALSE)
outputDir = 'outputs'
dir.create(outputDir)

#---------------------Main Loop--------------------------------------
for (scenName in scenDirs) {
  print(paste0("Working with ",scenName))
  tdata_raw = tibble(read.csv(paste0(scenName,'/analyzedTOF.csv')))
  nSub = length(tdata_raw$sim.id)
  
  # Who recovered sponatnously
  spntRec = tibble(
    nSpntRec = dim(tdata_raw |> filter(spntRecover==TRUE))[1],
    pSpntRec = nSpntRec / nSub)
  
  
  # exclude thse who didnt and actually never recovered
  tdata = tdata_raw |> 
    filter((spntRecover==FALSE) & (didRecover==TRUE))
  
  t90_q = quantileMe(tdata$t90,"t90") 
  
  RNB = tdata_raw |> filter((isRNB==TRUE) & (isDip==FALSE))
  npRNB = tibble(
    nRNB = dim(RNB)[1],
    pRNB = nRNB / dim(tdata)[1])
  
  initRNB_q = quantileMe(RNB$tIniRNB, "initRNB")
  durRNB_q = quantileMe(RNB$durRNB, "durRNB")
  lowestRNB_q = quantileMe(RNB$lowestTOF, "lowestRNB")
  
  DIP = tdata_raw |> filter((isDip==TRUE) & (spntRecover==FALSE))
  npDIP = tibble(
    nDIP = dim(DIP)[1],
    pDIP = nDIP / dim(tdata)[1] )
  
  initDIP_q = quantileMe(DIP$tDip, "initDIP")
  #lowestDIP_q = quantileMe(DIP$lowestTOF, "lowestDIP")
  
  scenLetter = substr(scenName,nchar(scenName),nchar(scenName))
  sumup = tibble(scen=scenLetter,
                     cbind(spntRec, t(t90_q),
                           npRNB, t(initRNB_q), t(durRNB_q), t(lowestRNB_q),
                           npDIP, t(initDIP_q)))
  if (!exists("scenSumup")) {
    scenSumup = sumup 
  }else {
    scenSumup = rbind(scenSumup, sumup)
  }
  
  #----------------PKPD params-----------------------------------
  #Build a table with percentiles 0,50,95 for parameter and scenario
  params = tibble(read.csv(paste0(scenName,'/params.csv'))) |>
    select(all_of(targetPars))
  # Sort events in 0(normal) 1(spntRec) 2(RNB) 3(dip) 4(neverRec)
  idEvents = tdata_raw |> 
    rowwise() |> 
    mutate(whichEv = whichEvent(spntRecover,isRNB,isDip,didRecover)) |>
    select(whichEv)
  
  # Add as new column to params
  params = params |>
    mutate(whichEv = idEvents$whichEv)
  pars_grp = params|> group_by(whichEv) |> 
    summarise(across(targetPars[2:9],
                     list(p5=~quantile(.,.05, na.rm=TRUE),
                          p50=~quantile(.,.5, na.rm=TRUE),
                          p95=~quantile(.,.95, na.rm=TRUE)
                          )))
  letterCol = tibble(scen=rep(scenLetter,dim(pars_grp)[1]))
  pars_grp = tibble(cbind(letterCol, pars_grp))
  
  if (!exists("parsSumup")) {
    parsSumup = pars_grp 
  }else {
    parsSumup = rbind(parsSumup, pars_grp)
  }

  
} #                        End of Main Loop

# Write sumup files------------------------------------------
write.csv(scenSumup, 	paste0(outputDir,'/TOFsumup.csv'), 	
          row.names = FALSE)
write.csv(parsSumup, 	paste0(outputDir,'/PARSsumup.csv'), 	
          row.names = FALSE)


# Create a human readable TOF data
TOFprettyPrinted = scenSumup |>
  rowwise() |>
  mutate(spntRecovered = npMe( nSpntRec,  pSpntRec),
         t90           = qMe( t90_p5, t90_p50, t90_p95),
         nRNB          = npMe( nRNB, pRNB),
         initRNB       = qMe( initRNB_p5, initRNB_p50, initRNB_p95),
         durRNB        = qMe( durRNB_p5, durRNB_p50, durRNB_p95),
         lowestRNB     = qMe( lowestRNB_p5, lowestRNB_p50, lowestRNB_p95),
         nDIP          = npMe( nDIP,  pDIP),
         initDIP      = qMe( initDIP_p5, initDIP_p50, initDIP_p95)) |>
  select(c(scen, spntRecovered,
           t90,
           nRNB,
           initRNB,
           durRNB,
           lowestRNB,
           nDIP,
           initDIP))

write.csv(TOFprettyPrinted, paste0(outputDir,'/TOFsumup_h.csv'), 
          row.names = FALSE)



