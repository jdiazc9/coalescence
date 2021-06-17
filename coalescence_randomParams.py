# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 09:57:19 2021

@author: juandc
"""



# %%
### INITIALIZATION

# reset variables etc.
from IPython import get_ipython
get_ipython().magic('reset -sf')

# libraries
import community_simulator
import community_simulator.usertools
import community_simulator.visualization
import numpy as np
#import pandas as pd
# import matplotlib
#import matplotlib.pyplot as plt
# import seaborn as sns
# colors = sns.color_palette()
# %matplotlib inline

import random
# import math
import copy
import os

# check if running on Windows
if os.name == 'nt':
    par = False
else:
    par = True

import time
start_time = time.time()

np.random.seed(1)
random.seed(1)



# %%
### USER-DEFINED FUNCTIONS

# stabilize communities through serial passaging
def stabilizeCommunities(plate):
    
    # parameters
    dilution_factor = 1/100
    n_transfers = 20
    time_per_cycle = 1
    # extinction_threshold = 1e-4 # we don't use this because we don't use StadyState()
    
    # run stabilization protocol
    Ntraj, Rtraj = plate.RunExperiment(f = np.diag(np.ones(plate.n_wells))*dilution_factor,
                                       T = time_per_cycle,
                                       npass = n_transfers,
                                       compress_resources=False,
                                       compress_species=True)
    
    return([plate.N,plate.R])

# Bray-Curtis similarity
# p and q are lists (of absolute abundances, sum need not be 1)
def bray_curtis(p,q):

    # C_pq: sum of lesser values of common species
    c = sum([min(p[i],q[i]) for i in range(len(p)) if p[i]>0 and q[i]>0])
    
    # Bray-Curtis similarity
    bc = 2*c/(sum(p) + sum(q))
    return bc

# similarity index
# p: invasive
# q: resident
# pq: coalesced
def mysim(p,q,pq):
    
    # if any of the vectors is empty, return nan (this could happen e.g. if maintenance costs are too high)
    if sum(p)==0 or sum(q)==0 or sum(pq)==0:
        return np.nan
    
    # otherwise...
    else:
    
        # normalization
        p = [p_i/sum(p) for p_i in p]
        q = [q_i/sum(q) for q_i in q]
        pq = [pq_i/sum(pq) for pq_i in pq]
        
        # relative similarity
        return bray_curtis(p,pq)/(bray_curtis(p,pq) + bray_curtis(q, pq))



# %%
### GENERATE RANDOM PARAMETERS, RUN MODEL, GET R^2, REPEAT

# list of R^2
r_squared = ['NA' for i in range(100)]

# repeat simulation 100 times
for n_sim in range(100):
    
    # print progress
    print('Running simulation %s out of 100' % (n_sim+1))
    time.sleep(0.1)
    
    # general assumptions (common to all simulations)
    assumptions = community_simulator.usertools.a_default.copy()
    assumptions['n_wells'] = 50 # number of communities
    
    assumptions['response'] = 'type I'
    assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
    assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
    assumptions['supply'] = 'off'
    assumptions['R0_food'] = 1000
    
    assumptions['c0'] = 0.0 # background consumption rate in binary model
    assumptions['c1'] = 1.0 # specific consumption rate in binary model
    
    assumptions['rs'] = 0.0 # control parameter for resource secretion: 0 means random secretions, 1 means species only secrete resources they can consume (relevant only if 'metabolism' is 'specific')
    
    # randomly sample variable parameters
    assumptions['S'] = random.randint(10,90) # number of species sampled at initialization -- uniformly distributed between 10 and 90
    
    n_families = random.randint(1,5) # number of species families -- uniformly distributed between 1 and 5
    n_species_per_family = random.randint(400,1200) # number of species per family -- uniformly distributed between 400 and 1200
    assumptions['SA'] = [n_species_per_family] * n_families
    assumptions['Sgen'] = random.randint(100,380) # number of generalists -- uniformly distributed between 100 and 380
    
    n_resource_classes = random.randint(3,8) # number of resource classes -- uniformly distributed between 1 and 5
    n_resources_per_class = random.randint(3,17) # number of resources per class -- uniformly distributed between 3 and 17
    assumptions['MA'] = [n_resources_per_class] * n_resource_classes
    
    #assumptions['l'] = [random.uniform(0.45,0.95) for i in range(sum(assumptions['MA']))] # leakage fraction -- uniformly distributed between 0.45 and 0.95, resource-specific
    assumptions['l'] = random.uniform(0.45,0.95) # leakage fraction -- uniformly distributed between 0.45 and 0.95, resource-unspecific
    
    #assumptions['m'] = [random.uniform(0,0.2) for i in range(sum(assumptions['SA']) + assumptions['Sgen'])] # species maintenance costs -- uniformly distributed between 0 and 0.2, species-specific
    assumptions['m'] = random.uniform(0,0.2) # species maintenance costs -- uniformly distributed between 0 and 0.2, species-unspecific
    
    assumptions['q'] = random.uniform(0.5,1) # preference strength (0 for generalists, 1 for specialists) -- uniformly distributed between 0.5 and 0.9
    assumptions['sigc'] = random.uniform(2,4) # standard deviation of sum of consumption rates for Gaussian and Gamma models -- uniformly distributed between 2 and 4
    assumptions['sparsity'] = random.uniform(0.05,0.95) # variability in secretion fluxes among resources (must be less than 1) -- uniformly distributed between 0.05 and 0.95
    assumptions['fs'] = random.uniform(0.05,0.45) # fraction of secretion flux to resources of the same type as the consumed one -- uniformly distributed between 0.05 and 0.45
    assumptions['fw'] = random.uniform(0.05,0.45) # fraction of secretion flux to waste resources -- uniformly distributed between 0.05 and 0.45
    assumptions['metabolism'] = ['common','specific'][random.randint(0,1)] # common or species-specific D matrix -- chosen arbitrarily with equal probabilities
    
    # parameters (make new matrices, etc)
    params = community_simulator.usertools.MakeParams(assumptions)
    
    # dynamics
    def dNdt(N,R,params):
        return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)
    def dRdt(N,R,params):
        return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)
    
    dynamics = [dNdt,dRdt]
    
    
    
    # %%
    ### SEED AND STABILIZE COMMUNITIES
    
    # make species pools
    pool_overlap = 0.0 # 0: no overlap, 1: full overlap (same pool)
    n_species = sum(assumptions['SA']) + assumptions['Sgen']
    n_pool = int(int(n_species/2) + pool_overlap*int(n_species/2))
    
    n = list(range(n_species))
    random.shuffle(n)
    
    pools = [n[:n_pool],
             n[-n_pool:]]
    pools[0].sort()
    pools[1].sort()
    
    # function to sample from pool
    def sampleFromPool(pool_number):
        N = np.zeros((n_species,assumptions['n_wells']))
        for i in range(assumptions['n_wells']):
            Ni = random.sample(pools[pool_number],assumptions['S'])
            N[Ni,i] = 1
        return N
            
    # generic initial state (will modify this later for specific plate settings)
    init_state = community_simulator.usertools.MakeInitialState(assumptions)
    
    # initialize resident plate
    N0_resident = init_state[0]
    N0_resident.loc[:,:] = sampleFromPool(0)
    init_state_resident = (N0_resident,
                           init_state[1])
    resident_plate = community_simulator.Community(init_state_resident,
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # initialize invasive plate
    N0_invasive = init_state[0]
    N0_invasive.loc[:,:] = sampleFromPool(1)
    init_state_invasive = (N0_invasive,
                           init_state[1])
    invasive_plate = community_simulator.Community(init_state_invasive,
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # stabilize plates
    N_resident, R_resident = stabilizeCommunities(resident_plate)
    N_invasive, R_invasive = stabilizeCommunities(invasive_plate)
    
    # make and stabilize coalescence plate
    N0_coalescence = (N_resident + N_invasive)/2*(1/100) # added dilution factor
    init_state_coalescence = (N0_coalescence,
                              init_state[1]) # resource init_state same as before
    coalescence_plate = community_simulator.Community(init_state_coalescence,
                                                      dynamics,
                                                      params,
                                                      parallel=par,
                                                      scale=1e6)
    N_coalescence, R_coalescence = stabilizeCommunities(coalescence_plate)
    
    
    
    # %%
    ### MONOCULTURES
    
    # find most abundant species in each resident/invasive community
    dominants_resident = list(list(zip(*N_resident.idxmax()))[1])
    dominants_invasive = list(list(zip(*N_invasive.idxmax()))[1])
    all_species = list(list(zip(*N_resident.index))[1])
    
    # initialize resident monoculture plate
    init_state_monoculture = community_simulator.usertools.MakeInitialState(assumptions)
    N0_monoculture_resident = init_state_monoculture[0]
    N0_monoculture_resident.iloc[:,:] = 0
    N0_monoculture_invasive = copy.deepcopy(N0_monoculture_resident)
    
    # set dominant species' abundance to 1 for each well
    for i in range(len(dominants_resident)):
        n = np.where(np.isin(all_species,dominants_resident[i]))[0][0]
        N0_monoculture_resident.iloc[n,i] = 1
        
    for i in range(len(dominants_invasive)):
        n = np.where(np.isin(all_species,dominants_invasive[i]))[0][0]
        N0_monoculture_invasive.iloc[n,i] = 1
    
    # make plates
    init_state_monoculture_resident = (N0_monoculture_resident,
                                       init_state_monoculture[1])
    monoculture_resident_plate = community_simulator.Community(init_state_monoculture_resident,
                                                               dynamics,
                                                               params,
                                                               parallel=par,
                                                               scale=1e6)
    
    init_state_monoculture_invasive = (N0_monoculture_invasive,
                                       init_state_monoculture[1])
    monoculture_invasive_plate = community_simulator.Community(init_state_monoculture_invasive,
                                                               dynamics,
                                                               params,
                                                               parallel=par,
                                                               scale=1e6)
    
    # stabilize monocultures (we just allow to grow for a single cycle)
    monoculture_resident_plate.Propagate(1,compress_resources=False,compress_species=True)
    monoculture_invasive_plate.Propagate(1,compress_resources=False,compress_species=True)
    
    
    
    # %%
    ### PAIRWISE COMPETITION
    
    # make plate (rescale so most abundant species at t=0 has abundance 1 (times the scaling factor of the plate))
    N0_pairwise = (monoculture_resident_plate.N + monoculture_invasive_plate.N)
    N0_pairwise = N0_pairwise/N0_pairwise.max()
    init_state_pairwise = (N0_pairwise,
                           init_state[1])
    pairwise_plate = community_simulator.Community(init_state_pairwise,
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # stabilize
    N_pairwise, R_pairwise = stabilizeCommunities(pairwise_plate)
    
    
    
    # %%
    ### STATISTICS
    
    # get relative similarity between invasive and coalesced communities
    q = ['NA' for i in range(assumptions['n_wells'])]
    
    for i in range(assumptions['n_wells']):
        
        # community compositions as lists
        R = N_resident.iloc[:,i].tolist()
        I = N_invasive.iloc[:,i].tolist()
        C = N_coalescence.iloc[:,i].tolist()
        
        q[i] = mysim(I,R,C)
        
    # get fraction in pairwise competition
    f_pairwise = ['NA' for i in range(assumptions['n_wells'])]
    f_pairwise_null = ['NA' for i in range(assumptions['n_wells'])]
    F = N_pairwise/N_pairwise.sum()
    F_null = N0_pairwise/N0_pairwise.sum()
    F = F.droplevel(0)
    F_null = F_null.droplevel(0)
    for i in range(assumptions['n_wells']):
        f_pairwise[i] = F.loc[dominants_invasive[i]][i]
        f_pairwise_null[i] = F_null.loc[dominants_invasive[i]][i]
    
    
    
    # %%
    ### FIT AND SAVE R^2
    
    r_squared[n_sim] = np.corrcoef(f_pairwise,q)[0,1]**2
    
   
    
    
    
    
    
    
    
    
    
    
    
