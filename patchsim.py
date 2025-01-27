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
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)

def read_config(config_file):
    config_df = pd.read_csv(config_file,delimiter='=',names=['key','val'])
    configs = dict(zip(config_df.key,config_df.val))
    if 'Model' not in configs.keys():
        configs['Model'] = 'Mobility'
    return configs


def load_patch(configs):
    patch_df = pd.read_csv(configs['PatchFile'],names=['id','pops'],
                            delimiter=' ',dtype={'id':str,'pops':int})
    patch_df.sort_values('id',inplace=True)

    logger.info('Loaded patch attributes')
    return patch_df

def load_params(configs,patch_df):
    params = {}
    params['T'] = int(configs['Duration'])
    try:
        #params['beta'] = np.repeat(float(configs['ExposureRate']),len(patch_df))
        params['beta'] = np.ones([len(patch_df),params['T']])*float(configs['ExposureRate'])
        params['alpha'] = float(configs['InfectionRate'])
        params['gamma'] = float(configs['RecoveryRate'])
        logger.info('Loaded disease parameters from Config')

    except:
        params['beta'] = np.zeros([len(patch_df),params['T']])
        params['alpha'] = np.repeat(0.0,len(patch_df))
        params['gamma'] = np.repeat(0.0,len(patch_df))
        logger.info('No parameter values in Config. Setting default to 0.')

    try:
        param_df = pd.read_csv(configs['ParamFile'], delimiter=' ',dtype={0:str},header=None).set_index(0).fillna(method='ffill',axis=1)
        patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
        param_df['Id_int'] = param_df.index.map(patch_idx)
        param_df.sort_values('Id_int',inplace=True)
        params['beta'][param_df.Id_int.values,:] = param_df[param_df.columns.difference(['Id_int'])].values

        logger.info('Loaded disease parameters from ParamFile')
    except:
        logger.info('No ParamFile loaded')
        pass

    try:
        params['vaxeff'] = float(configs['VaxEfficacy'])
    except:
        params['vaxeff'] = 1.0

    return params

def load_seed(configs,params,patch_df):
    try:
        seed_df = pd.read_csv(configs['SeedFile'],delimiter=' ',names=['Day','Id','Count'],dtype={'Id':str})
    except:
        empty_seed = np.ndarray((params['T'],len(patch_df)))
        empty_seed.fill(0.0)

        logger.info('Continuing without seeding')
        return empty_seed

    patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
    seed_df['Id_int'] = seed_df.Id.apply(lambda x: patch_idx[x])
    seed_df = seed_df.pivot(index='Day',columns='Id_int',values='Count').fillna(0)
    seed_df = seed_df.reindex(index=range(params['T']),columns = range(len(patch_df))).fillna(0)

    logger.info('Loaded seeding schedule')
    return seed_df.values

def load_vax(configs,params,patch_df):
    try:
        vax_df = pd.read_csv(configs['VaxFile'],delimiter=' ',
                    names=['Day','Id','Count'],dtype={'Id':str,'Count':int})
    except:
        empty_vax = np.ndarray((params['T'],len(patch_df)))
        empty_vax.fill(0.0)
        return empty_vax

    try:
        vax_delay = int(configs['VaxDelay'])
    except:
        vax_delay = 0

    patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
    vax_df['Id_int'] = vax_df.Id.apply(lambda x: patch_idx[x])
    vax_df['Delayed_Day'] = vax_df['Day'] + vax_delay

    vax_df = vax_df.pivot(index='Delayed_Day',columns='Id_int',values='Count').fillna(0)
    vax_df = vax_df.reindex(index=range(params['T']),columns = range(len(patch_df))).fillna(0)
    return vax_df.values.astype(int)

def load_Theta(configs, params, patch_df):
    theta_df = pd.read_csv(configs['NetworkFile'],names=['src_Id','dest_Id','theta_index','flow'],
                            delimiter=' ',dtype={'src_Id':str, 'dest_Id':str})

    if (configs['NetworkType']=='Static') & (len(theta_df.theta_index.unique())!=1):
        logger.info("Theta indices mismatch. Ensure NetworkType=Static.")
    if (configs['NetworkType']=='Weekly') & (len(theta_df.theta_index.unique())!=53):
        logger.info("Theta indices mismatch. Ensure NetworkType=Weekly.")
    if (configs['NetworkType']=='Monthly') & (len(theta_df.theta_index.unique())!=12):
        logger.info("Theta indices mismatch. Ensure NetworkType=Monthly.")


    patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
    try:
        theta_df['src_Id_int'] = theta_df.src_Id.apply(lambda x: patch_idx[x])
        theta_df['dest_Id_int'] = theta_df.dest_Id.apply(lambda x: patch_idx[x])
    except:
        logger.info("Ignoring flow entries for missing patches. Ensure all patches listed in PatchFile.")

    Theta_indices = theta_df.theta_index.unique()
    Theta = np.ndarray((len(Theta_indices),len(patch_df),len(patch_df)))

    for k in Theta_indices:
        theta_df_k = theta_df[theta_df.theta_index==k]
        theta_df_k = theta_df_k.pivot(index='src_Id_int',columns='dest_Id_int',values='flow').fillna(0)
        theta_df_k = theta_df_k.reindex(index=range(len(patch_df)),columns = range(len(patch_df))).fillna(0)
        Theta[int(k)] = theta_df_k.values

    logger.info('Loaded temporal travel matrix')
    return Theta

def patchsim_step(State_Array,patch_df,configs,params,theta,seeds,vaxs,t,stoch):
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
        ## vaccination for day t
        actual_vax = np.minimum(vaxs[t]*params['vaxeff'],S[t])
        S[t] =  S[t] - actual_vax
        V[t] = V[t] + actual_vax

        N = patch_df.pops.values

        ## Computing force of infection

        if configs['Model'] == 'Mobility':
            N_eff = theta.T.dot(N)
            I_eff = theta.T.dot(I[t])
            beta_j_eff = np.nan_to_num(np.multiply(np.divide(I_eff,N_eff),params['beta'][:,t]))
            inf_force = theta.dot(beta_j_eff)


        elif configs['Model'] == 'Force':
            beta_j_eff = np.nan_to_num(np.multiply(np.divide(I[t],N),params['beta'][:,t]))
            inf_force = theta.T.dot(beta_j_eff)

        ## New exposures during day t
        new_inf = np.multiply(inf_force,S[t])

        # print('=======Day {}========'.format(t))
        # print("theta.T")
        # print(theta.T)
        # print("I[t], N, beta",I[t],N,params['beta'][:,t])
        # print("beta_j_eff", beta_j_eff)
        # print("inf_force",inf_force)
        # print("S[t], new_inf", S[t],new_inf)


        S[t+1] = S[t] - new_inf
        E[t+1] = new_inf + np.multiply(1 - params['alpha'],E[t])
        I[t+1] = np.multiply(params['alpha'],E[t]) + np.multiply(1 - params['gamma'],I[t])
        R[t+1] = R[t] + np.multiply(params['gamma'],I[t])
        V[t+1] = V[t]


def epicurves_todf(configs,patch_df,State_Array):
    S,E,I,R,V = State_Array ## Aliases for the State Array

    try:
        scaling = float(configs['ScalingFactor'])
    except:
        scaling = 1

    rounding_bool = True
    try:
        if configs['OutputFormat'] == 'Whole':
            rounding_bool = True
        if configs['OutputFormat'] == 'Fractional':
            rounding_bool = False
    except:
        pass

    out_df = pd.DataFrame(columns = range(int(configs['Duration']) + 1))

    for i in range(len(patch_df)):
        net_sus = S[:,i]+V[:,i]
        if configs['LoadState']=='False':
            net_sus = np.lib.pad(net_sus,(1,0),'constant',constant_values=(patch_df.pops.values[i],))
        new_exposed = np.abs(np.diff(net_sus))


        epicurve = [int(x*scaling) if rounding_bool else (x*scaling) for x in new_exposed]
        out_df.loc[patch_df.id.values[i]] = epicurve

    return out_df


def write_epicurves(configs,patch_df,State_Array):
    S,E,I,R,V = State_Array ## Aliases for the State Array
    f = open(configs['OutputFile'],'w')

    try:
        scaling = float(configs['ScalingFactor'])
    except:
        scaling = 1

    rounding_bool = True
    try:
        if configs['OutputFormat'] == 'Whole':
            rounding_bool = True
        if configs['OutputFormat'] == 'Fractional':
            rounding_bool = False
    except:
        pass

    for i in range(len(patch_df)):
        net_sus = S[:,i]+V[:,i]
        if configs['LoadState']=='False':
            net_sus = np.lib.pad(net_sus,(1,0),'constant',constant_values=(patch_df.pops.values[i],))
        new_exposed = np.abs(np.diff(net_sus))

        if rounding_bool:
            epicurve = ' '.join([str(int(x*scaling)) for x in new_exposed])
        else:
            epicurve = ' '.join([str(x*scaling) for x in new_exposed])

        f.write('{} {}\n'.format(patch_df.id.values[i], epicurve))
    f.close()

def run_disease_simulation(configs,patch_df=None,params=None,Theta=None,seeds=None,vaxs=None,return_epi=False,write_epi=False):
    try:
        handler = logging.FileHandler(configs['LogFile'], mode='w')
        for hdlr in logger.handlers[:]:  # remove the existing file handlers
            if isinstance(hdlr,logger.FileHander):
                logger.removeHandler(hdlr)

        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
    except:
        handler = logging.NullHandler()
        logger.addHandler(handler)

    logger.info('Starting PatchSim')
    start = time.time()

    if configs['Model'] not in ['Mobility','Force']:
        logger.info('Invalid Model for PatchSim')
        logger.removeHandler(handler)
        return
    else:
        logger.info('Operating PatchSim under {} Model'.format(configs['Model']))

    if patch_df is None:
        patch_df = load_patch(configs)

    if params is None:
        params = load_params(configs,patch_df)

    if Theta is None:
        Theta = load_Theta(configs, params, patch_df)

    if seeds is None:
        seeds = load_seed(configs,params,patch_df)

    if vaxs is None:
        vaxs = load_vax(configs,params,patch_df)

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
    for t in range(params['T']):
        curr_date = ref + timedelta(days=t+int(configs['StartDate']))
        curr_week = int(curr_date.strftime("%U"))
        curr_month = int(curr_date.strftime("%m"))

        if configs['NetworkType']=='Static':
            patchsim_step(State_Array,patch_df,configs,params,Theta[0],seeds,vaxs,t,stoch)

        if configs['NetworkType']=='Weekly':
            patchsim_step(State_Array,patch_df,configs,params,Theta[curr_week-1],seeds,vaxs,t,stoch)

        if configs['NetworkType']=='Monthly':
            patchsim_step(State_Array,patch_df,configs,params,Theta[curr_month-1],seeds,vaxs,t,stoch)

    if configs['SaveState'] == 'True':
        logger.info('Saving StateArray to File')
        np.save(configs['SaveFile'],State_Array[:,-1,:])

    elapsed = time.time() - start
    logger.info('Simulation complete. Time elapsed: {} seconds.'.format(elapsed))
    logger.removeHandler(handler)

    if (write_epi==False)&(return_epi==False):
        return int(sum(R[-1,:]))
    else:
        if (write_epi==True):
            write_epicurves(configs,patch_df,State_Array)
            return

        if (return_epi==True):
            return epicurves_todf(configs,patch_df,State_Array)
