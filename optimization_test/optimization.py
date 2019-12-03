#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import copy
import logging
import simulation as sim
# from simulation import find_month,load_travel_network,patch_ODE_step
logger = logging.getLogger("__name__")

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


# Convert dict type vax_supply to np.array representation.
def vax_supply_to_array(supply_dict,T):
    supply = np.zeros(T + 1)
    for t in supply_dict.keys():
        if t <= T:
            supply[t] = supply_dict[t]
        else:
            logger.warning('supply on day {} comes after simulation ends'
                           .format(t))
    return supply


# Assign Q to given state: Q is divided among counties proportionally to
#  their population sizes
# vaxs: vax assignment schedule at county level
#  {day:{county_id:amount}}
# state: 2-digit state_id
# patch: county data
def assign_to_state_proportional(vaxs,state,patch,Q,day):
    assignment = copy.deepcopy(vaxs)
    patch_ids = patch.id
    patch_pops = patch.pops
    patch_count = len(patch_ids.values)
    state_popsize = 0.
    assignment[day] = []
    for i in range(patch_count):
        if patch_ids[i][0:2] == state:
            state_popsize += patch_pops[i]
    # now divide proportionally
    for i in range(patch_count):
        assignment[day].append(0.0)
        if patch_ids[i][0:2] == state:
            county = patch_ids[i]
            amount = int(round(patch_pops[i] / state_popsize * Q))
            if day not in assignment.keys():
                assignment[day] = {county:amount}
            else:
                assignment[day][i] += amount
                # if county not in assignment[day].keys():
                #     assignment[day][county] = amount;
                # else:
                #     assignment[day][county] += amount
    return assignment


# Check availability of vax for specific day.
# inventory[t]: amount added to inventory on day t
def check_avail(inventory,day,Q):
    return min(int(sum(inventory[0:day+1])),Q)

# Check feasibility of allocation by comparing
# against the state population caps

def check_popcaps(vaxs,state,patch,Q):
    patch_ids = patch.id
    patch_pops = patch.pops
    patch_count = len(patch_ids.values)
    state_popsize = 0
    for i in range(patch_count):
        if patch_ids[i][0:2] == state:
            state_popsize += patch_pops[i]
    state_alloc = 0
    for day in vaxs.keys():
        for i in range(patch_count):
            if patch_ids[i][0:2]==state:
                state_alloc+=vaxs[day][i]

    if (state_alloc+Q) > state_popsize:
        return False
    else:
        return True

# remove Q from most recent supply
def consume(supply,day,Q):
    to_consume = Q
    t = day
    while to_consume > 0 and t >= 0:
        if supply[t] > 0:
            consumption = min(supply[t],to_consume)
            supply[t] -= consumption
            to_consume -= consumption
        t -= 1
    if to_consume > 0:
        logger.error('supply is insufficient for consuming {} on day {}'.
                     format(Q,day))


# B is a dict {day:quantity_of_supply}
# L is step size
def greedy(cfg,B,L,temporal_granularity,look_ahead,f_vax_schedule):
    # import dispy
    # import os, socket
    # from sets import Set
    # host_set = Set()
    # node_file = open(os.environ['PBS_NODEFILE'],'r')
    # for line in node_file:
    #     host_set.add(line.strip())
    # hosts = list(host_set)
    # hosts.remove(socket.gethostname())
    #
    # cluster = dispy.JobCluster(sim.run_disease_simulation,nodes=hosts,
    #                            depends=[find_month,load_travel_network,
    #                                     patch_ODE_step])
    patch = sim.load_patch_attr(cfg)

    patch_ids = patch.id
    patch_pops = patch.pops
    patch_count = len(patch_ids.values)
    params = sim.load_disease_params(cfg, patch)
    seeds = sim.load_seed_schedule(cfg, params, patch)
    theta = sim.load_Theta(cfg, params, patch)
    T = params['T']
    budget = vax_supply_to_array(B,T)

    state_ids = sorted(list(set([id[0:2] for id in patch_ids.values])))
    #state_ids.remove('02') # Alaska
    #state_ids.remove('15') # Hawaii
    n_states = len(state_ids)
    u=1
    vaxs = {}
    configs = {"StartDate": 0, "NetworkType": "Weekly"}
    State_Array = None
    for day in range(1,T,temporal_granularity):
        Q = check_avail(budget,day,L)
        print("Day {}, Q {}".format(day, Q))
        ##Evaluating episize only until day+look_ahead duration
        if look_ahead!=-1:
            params['T'] = min(day+look_ahead,T)
        configs["StartDate"] = day
        state_arrs = [None] * n_states
        opt_j = 0
        while Q>0:
            print("    Budget Left {}".format(Q))
            jobs = []
            results = []
            logger.debug('Starting Greedy round {} on day {} looking ahead to day {}'.format(u,day,params['T']))
            logger.debug('{} units of vax available on day {}. Assigning...'.format(Q,day))

            for j in range(n_states):
                state = state_ids[j]
                ##Checking population limits
                if check_popcaps(vaxs,state,patch,Q):
                    x = assign_to_state_proportional(vaxs,state,patch,Q,day)
                    # job = cluster.submit(cfg,patch,params,seeds,x)
                    #job.id = (j,state)
                    # jobs.append(job)
                    episize, state_arrs[j] = sim.run_disease_simulation(cfg, patch_df=patch, params=params,Theta=theta, seeds=seeds, vaxs=x, input_state=State_Array, write_epi=)
                    results.append((j, state, episize))
                else:
                    logger.debug('Skipping state {} due to population limit'.format(state))

            # for job in jobs:
            #     j,state = job.id
            #     episize = job()
            #     results.append((j,state,episize))
            #     logger.debug('Vax {} to state FIPS {} on day {} with episize {}'.format(Q,state,day,episize))

            opt = sorted(results, key=lambda x: x[2])[0]
            opt_j,opt_state,opt_episize = opt

            logger.info('Round {} Best allocation: Vax {} to state FIPS {} on day {} with episize {}'.format(u,Q,opt_state,day,opt_episize))
            vaxs = assign_to_state_proportional(vaxs,state_ids[opt_j],patch,Q,day)
            consume(budget,day,Q)
            Q = check_avail(budget,day,L)
        u+=1
        State_Array = state_arrs[opt_j]
        logger.debug('Vax supply for day {} emptied. Continuing...'.format(day))
    logger.debug('All vaccines allocated!')
    sim.write_epicurves(cfg,patch,State_Array):

    for t in vaxs.keys():
        for i in range(len(patch_ids.values)):
            print('{} {} {}'.format(t,patch_ids.values[i],vaxs[t][i]))
            f_vax_schedule.write('%d %s %d\n' % (t,patch_ids.values[i],vaxs[t][i]))
