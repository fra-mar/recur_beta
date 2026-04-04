#!/usr/bin/env python
# coding: utf-8

# Kleijn model for interaction ROC-SUG 
# Based on Britsh Journal Clinical Pharmacology 2011; 72: 415-433
# Results have been compared to same model implemented in R library RxODE2
# https://www.omnicalculator.com/chemistry/molarity
# SUG molar mass: 2002.2 g/mol https://pubchem.ncbi.nlm.nih.gov/compound/sugammadex
# ROC molar mass: 529.8 g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Rocuronium

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
from kleijnModel_251104 import model
from paramsBuilder_251019 import buildParams
import gc

rng = np.random.default_rng(seed=1234)

args = sys.argv
if len(args) == 1:
    print('Script called without arguments for scenario number and nSubs')
else:
    # Which scenario (0...n scenarios)
    scenN = int(args[1])
    # N subjects
    nSub  = int(args[2])

# Load scenarios
df_scenarios = pd.read_csv('scenarios.csv')

# Follow up after SUG administration
follup = 3 # hours follow up after SUG administration

# Select the scenario numbered when the script was called
s = df_scenarios.iloc[scenN,:]
# 'simulations' folder should exist / created previously by job.sh script
scenDir = os.path.join('simulations',f'scenario_{df_scenarios.sceName[scenN]}')
os.mkdir(scenDir)

# Load variables for events (fe), and pkpd parameters (params)
params = buildParams(nSub, df_scenarios.durInf[scenN])


for i in range(0,nSub):
    # Unpack parameters for that individual
    ids, AGE, BW, RAC, SEV, infROC, durInf, d0ROC, catSUG = params.iloc[i,0:9]
    pars = params.iloc[i,9:]
    
    # Prepare a timeline
    # Two time spaces, one for ROC phase, other of SUG phase
    timeline  = np.linspace(0,durInf*60, int(durInf*60*6)) # for .1Hz srate
    timeline_ = np.linspace(0,follup*60, int(follup*60*6)) # for .1Hz srate

    # Both timelines concatenated
    TIMELINE  = np.hstack((timeline, timeline_ + timeline[-1]))

    # Total dROC in micromol (already incl BW), saved for later use
    dROC_tosave = d0ROC + infROC*60*durInf #x60 for µmol/min->µmol/h
    
    #-----------------------------ROC phase ------------------------------------
    
    # kind reminder, d0ROC, infROC Units are micromol
    
    # Initial values, ROC phase. 
    # AMOUNTS
    rocCentral, rocPerif = d0ROC, 0
    sugCentral, sugPerif = 0, 0
    # CONCENTRATIONS
    C1cx, C2cx           = 0, 0
    Ceff                 = 0

    
    # Initial values, ROC phase.
    inits = (rocCentral, rocPerif,
             sugCentral, sugPerif, 
             C1cx, C2cx,
             Ceff)  
    
    
    # Run the diff equations
    re = odeint(model, inits, timeline, args=(np.array(pars),infROC))
    ROCcentral, ROCperif, _, _, _, _, C_ROCeff = re.T # Capital letters to human-read them as np.arrays
    # Get last values as initial values for SUG phase
    next_inits = re[-1,:]

    # Calculated TOFr
    Emax = pars.E0
    tofR = pars.E0 - ( (Emax * C_ROCeff**pars.Hill) / (pars.EC50**pars.Hill + C_ROCeff**pars.Hill)  )
    tofR = np.nan_to_num(x=tofR, nan=pars.E0) #when high TOF line above can RuntimeWarning
        
    # -------------------------SUG Phase---------------------------------------
    lastTOFr = tofR[-1]
    
    if catSUG=='fixed':
        d0SUG = 200/ BW                        # mg/Kg
    elif catSUG == "adjusted":
        d0SUG = 2.0 if lastTOFr > 14 else 4.0  # mg/Kg
    else:
        print('found nothing about SUG dose')
    
    d0SUG = d0SUG*BW * .499451     # SUG dose microm, incl BW
    dSUG_tosave = d0SUG            # Spare for later comparison



    # Get last values will be inits for SUG step
    rocCentral, rocPerif, sugCentral, sugPerif, C1cx, C2cx, Ceff = next_inits 
    sugCentral = d0SUG  # Correct the sugCentral with the actual dose SUG 
    inits = (rocCentral, rocPerif,
             sugCentral, sugPerif,
             C1cx, C2cx, 
             Ceff) 
    
    re = odeint(model, inits, timeline_, args=(np.array(pars),0))

    ROCcentral_s, ROCperif_s, _, _, _, _, C_ROCeff_s = re.T

    # Calculate TOFr for SUG phase
    tofR_ = pars.E0 - ( (Emax * C_ROCeff_s**pars.Hill) /
                       (pars.EC50**pars.Hill +C_ROCeff_s**pars.Hill)  )
    
    tofR_ = np.nan_to_num(x=tofR_, nan=pars.E0) #when high TOF line above can RuntimeWarning
        
      # Transform amounts to concentrations

    if i == 0:
        C1ro = np.hstack(( ROCcentral, ROCcentral_s)).reshape(-1,1) /pars.V1ro 
        C2ro = np.hstack(( ROCperif,   ROCperif_s)).reshape(-1,1)   /pars.V2ro 
        Ceff_ro = np.hstack(( C_ROCeff,   C_ROCeff_s)).reshape(-1,1)  # Risk confusion with Ceff!!!
        TOFr = np.hstack((tofR, tofR_)).reshape(-1,1) 
        rSUG_remROC = d0SUG / (ROCcentral[-1] + ROCperif[-1])

        forRegression = np.array([[int(i+1), dROC_tosave,
                                   dSUG_tosave,
                                   ROCcentral[-1],
                                   ROCcentral[-1]/pars.V1ro,
                                   ROCperif[-1],
                                   ROCperif[-1]/pars.V2ro,
                                   rSUG_remROC]]
                                   )
    
    else:
        C1ro = np.hstack(   (C1ro, 
                             np.hstack(( ROCcentral, ROCcentral_s)).reshape(-1,1) /pars.V1ro)   )
        C2ro = np.hstack(   (C2ro,
                             np.hstack(( ROCperif,   ROCperif_s)).reshape(-1,1)   /pars.V2ro)   )
        Ceff_ro = np.hstack((Ceff_ro,
                             np.hstack(( C_ROCeff,   C_ROCeff_s)).reshape(-1,1) ))  # Risk confusion with Ceff!!!
        TOFr = np.hstack((TOFr,
                          np.hstack((tofR, tofR_)).reshape(-1,1))  )
        rSUG_remROC = d0SUG / (ROCcentral[-1] + ROCperif[-1])
        forRegression = np.vstack( (forRegression,
                                   np.array([[int(i+1), 
                                               dROC_tosave, 
                                               dSUG_tosave,
                                               ROCcentral[-1],
                                               ROCcentral[-1]/pars.V1ro,
                                               ROCperif[-1],
                                               ROCperif[-1]/pars.V2ro,
                                               rSUG_remROC]]) )
                                               )
          

  
    #------------------End of for loop#


ids = (np.arange(0,nSub) +1).astype('str')
ids = ['id'+x for x in ids]
ids = ",".join( ids )

np.savetxt(os.path.join(scenDir,'C1ro.csv'), 
           C1ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'C2ro.csv'), 
           C2ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'Ceff_ro.csv'), 
           Ceff_ro, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'TOFr.csv'), 
           TOFr, delimiter=',',header=ids, comments='')
np.savetxt(os.path.join(scenDir,'timeline.csv'), 
           TIMELINE, delimiter=',',header='time', comments='')
np.savetxt(os.path.join(scenDir,'forRegression.csv'),
           forRegression, delimiter=',', header='id,dROC,dSUG,A1ro,C1ro,A2ro,C2ro,rSUG_remROC', comments='')


# Save params

params_idx = np.array(params.index) +1
params.reset_index(drop=True, inplace=True)
params = pd.concat( (pd.DataFrame(data={'simId':params_idx}), params), axis=1)

params.to_csv(os.path.join(scenDir,'params.csv'), index=False)

print(f'\tDone with {scenDir}')
del C1ro, C2ro, Ceff_ro, TOFr, params, TIMELINE
gc.collect()


