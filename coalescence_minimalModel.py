# -*- coding: utf-8 -*-
"""
Created on Wed May 19 16:27:04 2021
@author: Juan
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
import pandas as pd
# import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette()
# %matplotlib inline

import random
import math
import copy
import os

# check if running on Windows
if os.name == 'nt':
    par = False
else:
    par = True

# set random and numpy's seed manually for reproducibility
#np.random.seed(0)
#random.seed(0)

import time
start_time = time.time()



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

# similarity index (relative Bray-Curtis similarity)
# p: invasive
# q: resident
# pq: coalesced
def mysim(p,q,pq):
    
    # normalization
    p = [p_i/sum(p) for p_i in p]
    q = [q_i/sum(q) for q_i in q]
    pq = [pq_i/sum(pq) for pq_i in pq]
    
    return bray_curtis(p,pq)/(bray_curtis(p,pq) + bray_curtis(q, pq))



# %%
### MODEL DEFINITION

# general assumptions
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 1 # number of communities
assumptions['S'] = 1
assumptions['SA'] = 4
assumptions['Sgen'] = 0 # 30 # number of generalists
assumptions['MA'] = 5 # [30, 30, 30] # number of resources per resource class
assumptions['l'] = 0.8 # leakage fraction

assumptions['response'] = 'type I'
assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)

assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species

# parameters
params = community_simulator.usertools.MakeParams(assumptions)

# dynamics
def dNdt(N,R,params):
    return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)

dynamics = [dNdt,dRdt]

# define competition experiment
#def competeMinimalCommunities(x, c_1inv_0 = 1, b_inv = 0, b_res = 0, c_2inv_alpha = 1, c_2res_gamma = 1):
    
x = 1.5
c_1res_0 = 1
b_inv = 0
b_res = 0
c_2inv_alpha = 1
c_2res_gamma = 1

# consumer matrix
params['c'].iloc[:,:] = np.array([[x, 0, b_inv, 0, 0],
                               [0, c_2inv_alpha, 0, 0, 0],
                               [c_1res_0, 0, 0, 0, b_res],
                               [0, 0, 0, c_2res_gamma, 0]])

# metabolic matrices
for i in range(len(params['D'])):
    params['D'][i].iloc[:,:] = np.identity(int(assumptions['MA']))

params['D'][0].iloc[:,0] = np.array([0, 1, 0, 0, 0])
params['D'][1].iloc[:,1] = np.array([0, 0, 1, 0, 0])
params['D'][2].iloc[:,0] = np.array([0, 0, 0, 1, 0])
params['D'][3].iloc[:,3] = np.array([0, 0, 0, 0, 1])



# %%
### COALESCENCE

# generic initial state
init_state = community_simulator.usertools.MakeInitialState(assumptions)

# invasive community
N0_inv = init_state[0].copy()
N0_inv.iloc[:,0] = 0
N0_inv.iloc[0:2,0] = 1
plate_inv = community_simulator.Community((N0_inv,init_state[1]),
                                          dynamics,
                                          params,
                                          parallel=par,
                                          scale=1e6)

# resident community
N0_res = init_state[0].copy()
N0_res.iloc[:,0] = 0
N0_res.iloc[2:5,0] = 1
plate_res = community_simulator.Community((N0_res,init_state[1]),
                                          dynamics,
                                          params,
                                          parallel=par,
                                          scale=1e6)

# stabilize primary communities
N_inv, R_inv = stabilizeCommunities(plate_inv)
N_res, R_res = stabilizeCommunities(plate_res)

# coalesced community
plate_coa = community_simulator.Community(((N_inv+N_res)/2,init_state[1]),
                                          dynamics,
                                          params,
                                          parallel=par,
                                          scale=1e6)

# stabilize coalesced community
N_coa, R_coa = stabilizeCommunities(plate_coa)





