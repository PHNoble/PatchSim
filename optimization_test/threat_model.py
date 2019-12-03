#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' PatchSim v1.2
Created and maintained by: Srini (srini@virginia.edu)
Date last modified: 6 Aug 2019
'''
import numpy as np
import pandas as pd
import logging
import time
import support_modules as sm
from datetime import datetime, timedelta

#Numerically stable softmax
def softmax(x):
    shiftx = x - np.max(x)
    exps = np.exp(shiftx)
    return exps / np.sum(exps)

def patchsim_step(State_Array,supply,vax_alloc,patch_df,params,theta,seeds,vaxs,t,stoch,logger):
    S,E,I,R,V = State_Array ## Aliases for the State Array

    ## seeding for day t (seeding implies S->I)
    actual_seed = np.minimum(seeds[t],S[t])
    S[t] = S[t] - actual_seed
    I[t] = I[t] + actual_seed

    if stoch:
        ## vaccination for day t
        max_SV = np.minimum(vaxs[t],S[t])
        actual_SV = np.random.binomial(max_SV.astype(int),params['vaxeff'])
        S[t] = S[t] - actual_SV
        V[t] = V[t] + actual_SV

        ## Computing force of infection
        N = patch_df.pops.values
        S_edge = np.concatenate([np.random.multinomial(S[t][x],theta[x]).reshape(1,len(N)) for x in range(len(N))],axis=0)
        E_edge = np.concatenate([np.random.multinomial(E[t][x],theta[x]).reshape(1,len(N)) for x in range(len(N))],axis=0)
        I_edge = np.concatenate([np.random.multinomial(I[t][x],theta[x]).reshape(1,len(N)) for x in range(len(N))],axis=0)
        R_edge = np.concatenate([np.random.multinomial(R[t][x],theta[x]).reshape(1,len(N)) for x in range(len(N))],axis=0)
        V_edge = np.concatenate([np.random.multinomial(V[t][x],theta[x]).reshape(1,len(N)) for x in range(len(N))],axis=0)
        N_edge = S_edge + E_edge + I_edge + R_edge + V_edge

        N_eff = N_edge.sum(axis=0)
        I_eff = I_edge.sum(axis=0)
        beta_j_eff = np.nan_to_num(params['beta'][:,t]*(I_eff/N_eff))

        actual_SE = np.concatenate([np.random.binomial(S_edge[:,x],beta_j_eff[x]).reshape(len(N),1) for x in range(len(N))],axis=1).sum(axis=1)
        actual_EI = np.random.binomial(E[t],params['alpha'])
        actual_IR = np.random.binomial(I[t],params['gamma'])

        S[t+1] = S[t] - actual_SE
        E[t+1] = E[t] + actual_SE - actual_EI
        I[t+1] = I[t] + actual_EI - actual_IR
        R[t+1] = R[t] + actual_IR
        V[t+1] = V[t]

    else:

        ## Computing force of infection
        N = patch_df.pops.values
        N_eff = theta.T.dot(N)
        I_eff = theta.T.dot(I[t])
        beta_j_eff = np.nan_to_num(np.multiply(np.divide(I_eff,N_eff),params['beta'][:,t]))
        inf_force = theta.dot(beta_j_eff)
        ## New exposures during day t
        new_inf = np.multiply(inf_force,S[t])

        ## vaccination for day t
        #Get the areas experiencing the highest force of infection, allocate them more
        #Vaccines based on predicted number of new exposures
        dis = softmax(new_inf)
        #Allocate vaxines
        supt = supply.get(t,0)
        actual_vax = vax_alloc[t]
        for i in range(len(patch_df)):
            actual_vax[i] = np.round(min(dis[i]*supt, S[t][i]))
        #move any leftover vaccines to the next day
        if supply.get(t+1, 0) == 0:
            supply[t+1] = supply.get(t, 0) - sum(actual_vax)
        else:
            supply[t+1] += supply.get(t, 0) - sum(actual_vax)

        S[t] =  S[t] - actual_vax
        V[t] = V[t] + actual_vax
        logger.debug('At time {}, allocated {} vaccines out of {}'.format(t, sum(actual_vax), supply.get(t,0)))

        ## Computing force of infection
        N = patch_df.pops.values
        N_eff = theta.T.dot(N)
        I_eff = theta.T.dot(I[t])
        beta_j_eff = np.nan_to_num(np.multiply(np.divide(I_eff,N_eff),params['beta'][:,t]))
        inf_force = theta.dot(beta_j_eff)

        S[t+1] = S[t] - new_inf
        E[t+1] = new_inf + np.multiply(1 - params['alpha'],E[t])
        I[t+1] = np.multiply(params['alpha'],E[t]) + np.multiply(1 - params['gamma'],I[t])
        R[t+1] = R[t] + np.multiply(params['gamma'],I[t])
        V[t+1] = V[t]




def run_disease_simulation(configs,supply,logger,patch_df=None,params=None,Theta=None,seeds=None,vaxs=None,return_epi=False,write_epi=False):

    logger.info('Starting PatchSim')
    start = time.time()

    if patch_df is None:
        patch_df = sm.load_patch(configs, logger)

    if params is None:
        params = sm.load_params(configs,patch_df, logger)

    if Theta is None:
        Theta = sm.load_Theta(configs, params, patch_df, logger)

    if seeds is None:
        seeds = sm.load_seed(configs,params,patch_df, logger)

    if vaxs is None:
        vaxs = sm.load_vax(configs,params,patch_df)

    logger.info('Initializing simulation run...')

    if 'RandomSeed' in configs.keys():
        np.random.seed(int(configs['RandomSeed']))
        stoch = True
        logger.info('Found RandomSeed. Running in stochastic mode...')
    else:
        stoch = False
        logger.info('No RandomSeed found. Running in deterministic mode...')

    dim = 5 ##Number of states (SEIRV)
    if stoch:
        State_Array = np.ndarray((dim,params['T']+1,len(patch_df))).astype(int)
    else:
        State_Array = np.ndarray((dim,params['T']+1,len(patch_df)))

    State_Array.fill(0)
    S,E,I,R,V = State_Array ## Aliases for the State Array

    if configs['LoadState'] =='True':
        State_Array[:,0,:] = np.load(configs['LoadFile'])
    else:
        S[0,:] = patch_df.pops.values


    ref = datetime.strptime('Jan 1 2017', '%b %d %Y') ##is a Sunday
    vax_alloc = np.ndarray((params['T']+1,len(patch_df)))
    for t in range(params['T']):
        curr_date = ref + timedelta(days=t+int(configs['StartDate']))
        curr_week = int(curr_date.strftime("%U"))
        curr_month = int(curr_date.strftime("%m"))

        if configs['NetworkType']=='Static':
            patchsim_step(State_Array,supply,vax_alloc,patch_df,params,Theta[0],seeds,vaxs,t,stoch, logger)

        if configs['NetworkType']=='Weekly':
            patchsim_step(State_Array,supply,vax_alloc,patch_df,params,Theta[curr_week-1],seeds,vaxs,t,stoch,logger)

        if configs['NetworkType']=='Monthly':
            patchsim_step(State_Array,supply,vax_alloc,patch_df,params,Theta[curr_month-1],seeds,vaxs,t,stoch,logger)

    if configs['SaveState'] == 'True':
        logger.info('Saving StateArray to File')
        np.save(configs['SaveFile'],State_Array[:,-1,:])

    elapsed = time.time() - start
    logger.info('Simulation complete. Time elapsed: {} seconds.'.format(elapsed))
    logger.debug('Vax allocation \n{}'.format(vax_alloc))
    if (write_epi==False)&(return_epi==False):
        return int(sum(R[-1,:]))
    else:
        if (write_epi==True):
            sm.write_epicurves(configs,patch_df,State_Array)
            f = open(configs['VaxSchedule'], 'w')
            for t in range(params['T']):
                for i in range(len(patch_df)):
                    if vax_alloc[t][i] > 0:
                        f.write('{} {} {}\n'.format(t, patch_df.id.values[i], vax_alloc[t][i]))
            f.close()
        if (return_epi==True):
            return sm.epicurves_todf(configs,patch_df,State_Array)
