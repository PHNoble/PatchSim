import numpy as np
import pandas as pd
import logging


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


# Read vax supply schedule from file.
# Each line of file should be:
#  day,supply_quantity
#
def load_vax_supply(supply_schedule_file):
    supply = {}
    f = open(supply_schedule_file,'r')
    for line in f:
        day,amount = line.strip().split(',')
        day,amount = int(day),int(amount)
        if day not in supply.keys():
            supply[day] = amount
        else:
            supply[day] += amount
    f.close()
    return supply

def read_config(config_file):
    config_df = pd.read_csv(config_file,delimiter='=',names=['key','val'])
    configs = dict(zip(config_df.key,config_df.val))
    return configs

def load_patch(configs, logger=None):
    patch_df = pd.read_csv(configs['PatchFile'],names=['id','pops'],
                            delimiter=' ',dtype={'id':str,'pops':int})
    patch_df.sort_values('id',inplace=True)
    if logger is not None:
        logger.info('Loaded patch attributes')
    return patch_df

def load_params(configs,patch_df,logger=None):
    params = {}
    params['T'] = int(configs['Duration'])
    try:
        #params['beta'] = np.repeat(float(configs['ExposureRate']),len(patch_df))
        params['beta'] = np.ones([len(patch_df),params['T']])*float(configs['ExposureRate'])
        params['alpha'] = float(configs['InfectionRate'])
        params['gamma'] = float(configs['RecoveryRate'])
        if logger is not None:
            logger.info('Loaded disease parameters from Config')

    except:
        params['beta'] = np.zeros([len(patch_df),params['T']])
        params['alpha'] = np.repeat(0.0,len(patch_df))
        params['gamma'] = np.repeat(0.0,len(patch_df))
        if logger is not None:
            logger.info('No parameter values in Config. Setting default to 0.')

    try:
        param_df = pd.read_csv(configs['ParamFile'], delimiter=' ',dtype={0:str},header=None).set_index(0).fillna(method='ffill',axis=1)
        patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
        param_df['Id_int'] = param_df.index.map(patch_idx)
        param_df.sort_values('Id_int',inplace=True)
        params['beta'][param_df.Id_int.values,:] = param_df[param_df.columns.difference(['Id_int'])].values
        if logger is not None:
            logger.info('Loaded disease parameters from ParamFile')
    except:
        if logger is not None:
            logger.info('No ParamFile loaded')
        pass

    try:
        params['vaxeff'] = float(configs['VaxEfficacy'])
    except:
        params['vaxeff'] = 1.0

    return params

def load_seed(configs,params,patch_df, logger=None):
    try:
        seed_df = pd.read_csv(configs['SeedFile'],delimiter=' ',names=['Day','Id','Count'],dtype={'Id':str})
    except:
        empty_seed = np.ndarray((params['T'],len(patch_df)))
        empty_seed.fill(0.0)
        if logger is not None:
            logger.info('Continuing without seeding')
        return empty_seed

    patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
    seed_df['Id_int'] = seed_df.Id.apply(lambda x: patch_idx[x])
    seed_df = seed_df.pivot(index='Day',columns='Id_int',values='Count').fillna(0)
    seed_df = seed_df.reindex(index=range(params['T']),columns = range(len(patch_df))).fillna(0)
    if logger is not None:
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

def load_Theta(configs, params, patch_df, logger=None):
    theta_df = pd.read_csv(configs['NetworkFile'],names=['src_Id','dest_Id','theta_index','flow'],
                            delimiter=' ',dtype={'src_Id':str, 'dest_Id':str})

    if (configs['NetworkType']=='Static') & (len(theta_df.theta_index.unique())!=1) & (logger != None):
        logger.info("Theta indices mismatch. Ensure NetworkType=Static.")
    if (configs['NetworkType']=='Weekly') & (len(theta_df.theta_index.unique())!=53) & (logger != None):
        logger.info("Theta indices mismatch. Ensure NetworkType=Weekly.")
    if (configs['NetworkType']=='Monthly') & (len(theta_df.theta_index.unique())!=12) & (logger != None):
        logger.info("Theta indices mismatch. Ensure NetworkType=Monthly.")


    patch_idx = dict(zip(patch_df.id.values,range(len(patch_df))))
    try:
        theta_df['src_Id_int'] = theta_df.src_Id.apply(lambda x: patch_idx[x])
        theta_df['dest_Id_int'] = theta_df.dest_Id.apply(lambda x: patch_idx[x])
    except:
        if logger is not None:
            logger.info("Ignoring flow entries for missing patches. Ensure all patches listed in PatchFile.")
        pass

    Theta_indices = theta_df.theta_index.unique()
    Theta = np.ndarray((len(Theta_indices),len(patch_df),len(patch_df)))

    for k in Theta_indices:
        theta_df_k = theta_df[theta_df.theta_index==k]
        theta_df_k = theta_df_k.pivot(index='src_Id_int',columns='dest_Id_int',values='flow').fillna(0)
        theta_df_k = theta_df_k.reindex(index=range(len(patch_df)),columns = range(len(patch_df))).fillna(0)
        Theta[int(k)] = theta_df_k.values
    if logger is not None:
        logger.info('Loaded temporal travel matrix')
    return Theta
