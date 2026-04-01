#!/usr/bin/env python
# coding: utf-8
# Build params for Kleijn model implementation in Python

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import shutil
import json
import sys

roc2microM = 1.887505
sug2microM = 0.499451
# Fraction of parent RocBromide which is rocuronium
rocBromide2roc = 1 #529.8/609.7 

# for reproducibility
rng = np.random.default_rng(seed=1234)

# For same format as R
def rnorm(mean, sd):
    return np.random.normal(mean,sd)


def buildParams(nSub, durInf):
    
    for i in range(nSub):
      AGE    = np.random.normal(69,10)
      CR     = 80 if AGE > 60 else 100
      BW     = np.random.normal(79,10)
      RAC   = np.random.randint(0,3)
      RAC   = 'Asian' if RAC==0 else 'nonAsian' # 33% will be Asian
      #SEV    = np.random.randint(0,2,1).astype('bool')[0] #Sevo True/False
      SEV = True
      #i_ROC  = np.random.randint(6,19)*5/100 #for i_ROC 0.3-0.9 by .05 mgKgH
      i_ROC  = np.random.randint(3,10)/10 #for i_ROC 0.3-0.9 by .1 mgKgH
      d0ROC  = 0.6  * BW            #induction dose ROC in mg (mgKg-1 *Kg)
      d0ROC  = d0ROC * rocBromide2roc * roc2microM #conversion to micromoles
      infROC = i_ROC * BW / 60                           # ROC mg per minute
      infROC = infROC * rocBromide2roc * roc2microM #conversion to micromoles
      # 20% cases will get SUG fixed, i.e. 200mg (1 vial) / BW
      catSUG = 'fixed' if np.random.random() <= .2 else 'adjusted'
      #d0SUG  = bolusSUG * BW            # Reversal dose SUG in mg (mgKg-1 *Kg)
      #d0SUG  = d0SUG * sug2microM  # Conversion mg to micromoles
      #tSUG   = durInf                # timestamp when SUG administered (hours)
      '''
      forEvents = pd.DataFrame(index=[0],
                               data={'d0ROC':d0ROC, 
                                     'infROC':infROC, 
                                     'durInf':durInf, 
                                     'd0SUG':d0SUG, 
                                     'tSUG':tSUG})
      '''
      CLage   = 1 + -0.00678 * (AGE - 43)
      tvCLro  = CLage * 0.269 *(BW/70)**0.75
      V1cr    = np.exp(-0.00143 * (CR-119))
      tvV1ro  = V1cr * 4.73 * (BW/70)**1
      Q2roRac = 1 if RAC=='nonAsian' else 1-0.212
      Q2ro    = Q2roRac * 0.279 * (BW/70)**0.75
      V2roAge = np.exp(0.00613 * (AGE-43))
      tvV2ro  = V2roAge * 6.76 * (BW/70)**1
      
      # for SUG
      CLsuBW  = 1 + 0.00378 * (BW-74.5)
      REN     = (2 * CR/(CR+119))**1.29
      tvCLsu  = CLsuBW * (REN * 0.093)
      V1suBW  = 1 + -0.00354 *(BW- 74.5)
      V1suRAC = 1 if RAC=='nonAsian' else 1-0.16
      tvV1su  = V1suBW * V1suRAC * 4.42 * (BW/70)**1
      tvQ2su  = 0.206 * (BW/70)**0.75
      V2suCR  = np.exp( -0.00305 * (CR-119))
      tvV2su  = V2suCR * 6.35 * (BW/70)**1
      
      # PD model ROC
      ke0Sev  = 1 if SEV==False else  1-0.567 # 1 no sevo,  1 +  -0.567 if sevo
      tvke0   = ke0Sev * 0.134 * (BW/70)**-0.25
      EC50sev = 1 if SEV==False else 1-0.395 # 1 no sevo,  1 +  -0.395 if sevo
      tvEC50  = EC50sev * 1.62 # microM
      #PD model SUG
      tvks    = np.exp(-3.43) * (BW/70)**-0.25
      tvHill  = 7.52
      tvE0    = 100 #1.04 * 100 #Arbitrary set to 100
      
      # eta
      eta_E0   = 0.01224571 # CV 11.1%
      eta_V1ro = 0.05509785 # CV% 23.8%
      eta_CLro = 0.09982362 # CV% 32.4%
      eta_V2ro = 0.09865367 # CV% 32.2%
      eta_ke0  = 0.16032217 #CV % 41.7%
      eta_EC50 = 0.06015486 # CV 24.9%
      eta_Hill = 0.15608110 # CV 41.1%
      eta_CLsu = 0.04895777 # CV 22.4%
      eta_ks   = 0.83273519 # CV 114%
      
      # Log-normal random distributed parameters
      nSub = nSub - 1
      # To understand np.sqrt(eta_xxx) look at https://nmhelp.tingjieguo.com/V/4 4.1.1 
      parsRnd = pd.DataFrame(data = {
        'id'  : i+1,
        'AGE'  : AGE,
        'BW'   : BW,
        'RAC' : RAC,
        'SEV' : SEV,
        'iROC' : i_ROC,
        'durInf' : durInf,
        'd0ROC' : d0ROC,
        'catSUG' : catSUG,
        #'d0SUG' : d0SUG,
        #'tSUG' : tSUG,
        'CLro' : tvCLro * np.exp(rnorm(mean=0, sd= np.sqrt(eta_CLro))),
        'Q2ro' : Q2ro,
        'V1ro' : tvV1ro * np.exp(rnorm(mean=0, sd= np.sqrt(eta_V1ro))),
        'V2ro' : tvV2ro * np.exp(rnorm(mean=0, sd= np.sqrt(eta_V2ro))),
        'V1su' : tvV1su,
        'V2su' : tvV2su,
        'CLsu' : tvCLsu * np.exp(rnorm(mean=0, sd= np.sqrt(eta_CLsu))),
        'Q2su' : tvQ2su,
        'ke0' : tvke0 * np.exp(rnorm(mean=0, sd= np.sqrt(eta_ke0))),
        'EC50' : tvEC50 * np.exp(rnorm(mean=0, sd= np.sqrt(eta_EC50))),
        'Hill' : tvHill * np.exp(rnorm(mean=0, sd= np.sqrt(eta_Hill))),
        'E0' : tvE0, # * np.exp(rnorm(nSub,mean=0, sd= np.sqrt(eta_E0))),
        'ks' : tvks * np.exp(rnorm(mean=0, sd= np.sqrt(eta_ks)))}, index=[0]
      )
      
      if i==0 :
          params = parsRnd.copy()
      else:
          params = pd.concat((params, parsRnd))
    params.reset_index(inplace=True, drop=True)
    return params


if __name__ == "__main__":
    params = buildParams(1000)
    print(params.head(50))
    print(params.describe().T)

