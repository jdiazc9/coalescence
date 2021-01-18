# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 15:24:54 2021

@author: Juan
"""



### INITIALIZIATION

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
# import math
# import copy

# set random and numpy's seed manually for reproducibility
np.random.seed(0)
random.seed(0)

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

# stack plot with bars
def myStackPlot(N):

    # check number of wells and species
    n_species = N.shape[0]
    n_wells = N.shape[1]
    
    # normalize each column to get fractions
    N = N/N.sum()
    
    # numpy data frame to list of lists (list of as many lists as species in total)
    N = N.values.tolist()
    
    # colormap (adjusted to number of species)
    random.seed(0)
    cmap = plt.get_cmap("gist_stern")
    cmap = cmap(np.linspace(0,1,num=n_species)).tolist()
    random.shuffle(cmap)
    
    alpha = 0.75
    cmap = [[227/255,0,0,alpha],
            [255/255,161/255,77/255,alpha],
            [51/255,153/255,237/255,alpha],
            [168/255,255/255,140/255,alpha]]
    
    # position in the x axis
    xpos = np.arange(n_wells).tolist()
    xticklabs = np.add(xpos,[1]).tolist()
    
    # define cumulative sum of abundances
    cumsum = np.zeros(n_wells).tolist()
    
    # initialize figure
    fig, ax = plt.subplots(1)
    
    # plot bars - iterate across species and stack
    for i in range(n_species):
        ax.bar(xpos,N[i],color=cmap[i],bottom=cumsum,width=0.5)
        cumsum = np.add(N[i],cumsum).tolist()
        
    ax.set_xticks(xpos)
    ax.set_xticklabels(xticklabs)
    ax.set_ylabel('Fraction')
    ax.set_ylim(0,1)
        
    return fig, ax



### MODEL DEFINITION

# general assumptions
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 1 # number of communities
assumptions['SA'] = [4] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 0 # number of generalists
assumptions['l'] = 0.6 # leakage fraction
assumptions['MA'] = [4] # [30, 30, 30] # number of resources per resource class

assumptions['response'] = 'type I'
assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)

assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['crossfeeding_strength'] = 1 # controls the secretions of the secondary invasive species into the resource that the dominant invasive can consume

# initialize upatke rate matrix
M = np.sum(assumptions['MA'])
T = len(assumptions['MA'])
S = np.sum(assumptions['SA'])+assumptions['Sgen']
F = len(assumptions['SA'])
resource_names = ['R'+str(k) for k in range(M)]
type_names = ['T'+str(k) for k in range(T)]
family_names = ['F'+str(k) for k in range(F)]
consumer_names = ['S'+str(k) for k in range(S)]
resource_index = [[type_names[m] for m in range(T) for k in range(assumptions['MA'][m])],
                  resource_names]
consumer_index = [[family_names[m] for m in range(F) for k in range(assumptions['SA'][m])]
                  +['GEN' for k in range(assumptions['Sgen'])],consumer_names]
c = pd.DataFrame(np.zeros((S,M)),columns=resource_index,index=consumer_index)

# custom uptake rates
c.loc[('F0','S0')] = np.array([1,0,0,0])
c.loc[('F0','S1')] = np.array([0,0,0.1,0])
c.loc[('F0','S2')] = np.array([0.8,1,0,0])
c.loc[('F0','S3')] = np.array([0,0,1,0])

# initialize metabolic matrices
Dlist = ['NA' for i in range(M)]

# dominant resident (S0)
D = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
D.loc[('T0','R0')] = np.array([0,0,1,0])
D.loc[('T0','R1')] = np.array([0,0,0,1])
D.loc[('T0','R2')] = np.array([0,0,0,1])
D.loc[('T0','R3')] = np.array([0,0,0,1])
Dlist[0] = D.T

# secondary resident (S1)
D = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
D.loc[('T0','R0')] = np.array([0,0,0,1])
D.loc[('T0','R1')] = np.array([0,0,0,1])
D.loc[('T0','R2')] = np.array([0,0,0,1])
D.loc[('T0','R3')] = np.array([0,0,0,1])
Dlist[1] = D.T

# dominant invasive (S2)
D = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
D.loc[('T0','R0')] = np.array([0,0,1,0])
D.loc[('T0','R1')] = np.array([0,0,0,1])
D.loc[('T0','R2')] = np.array([0,0,0,1])
D.loc[('T0','R3')] = np.array([0,0,0,1])
Dlist[2] = D.T

# secondary invasive (S3)
D = pd.DataFrame(np.zeros((M,M)),index=c.keys(),columns=c.keys())
D.loc[('T0','R0')] = np.array([0,0,0,1])
D.loc[('T0','R1')] = np.array([0,0,0,1])
D.loc[('T0','R2')] = np.array([0,assumptions['crossfeeding_strength'],0,1-assumptions['crossfeeding_strength']])
D.loc[('T0','R3')] = np.array([0,0,0,1])
Dlist[3] = D.T

D = Dlist

# parameters
assumptions['sampling'] = 'Gamma'
assumptions['S'] = 1
params = community_simulator.usertools.MakeParams(assumptions)
params['c'] = c
params['D'] = D

# dynamics
def dNdt(N,R,params):
    return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)

dynamics = [dNdt,dRdt]



### SEED AND STABILIZE PLATES

# generic initial state (will modify this later for specific plate settings)
init_state = community_simulator.usertools.MakeInitialState(assumptions)

# resident plate
N0_resident = init_state[0]
N0_resident.loc[:,:] = np.array([[1,1,0,0]]).T
init_state_resident = (N0_resident,
                       init_state[1])
resident_plate = community_simulator.Community(init_state_resident,
                                               dynamics,
                                               params,
                                               parallel=False,
                                               scale=1e6)
N_resident, R_resident = stabilizeCommunities(resident_plate)

# invasive plate
N0_invasive = init_state[0]
N0_invasive.loc[:,:] = np.array([[0,0,1,1]]).T
init_state_invasive = (N0_invasive,
                       init_state[1])
invasive_plate = community_simulator.Community(init_state_invasive,
                                               dynamics,
                                               params,
                                               parallel=False,
                                               scale=1e6)
N_invasive, R_invasive = stabilizeCommunities(invasive_plate)



# monoculture of dominant invasive
N0_monoculture = init_state[0]
N0_monoculture.loc[:,:] = np.array([[0,0,1,0]]).T
init_state_monoculture = (N0_monoculture,
                          init_state[1])
monoculture_plate = community_simulator.Community(init_state_monoculture,
                                               dynamics,
                                               params,
                                               parallel=False,
                                               scale=1e6)
N_monoculture, R_monoculture = stabilizeCommunities(monoculture_plate)

# dominant invasive species invading alone
N0_singleinv = init_state[0]
N0_singleinv.loc[:,:] = (N_resident + N_monoculture)/2
init_state_singleinv = (N0_singleinv,
                        init_state[1])
singleinv_plate = community_simulator.Community(init_state_singleinv,
                                               dynamics,
                                               params,
                                               parallel=False,
                                               scale=1e6)
N_singleinv, R_singleinv = stabilizeCommunities(singleinv_plate)

# full coalescence
N0_coalescence = init_state[0]
N0_coalescence.loc[:,:] = (N_resident + N_invasive)/2
init_state_coalescence = (N0_coalescence,
                        init_state[1])
coalescence_plate = community_simulator.Community(init_state_coalescence,
                                               dynamics,
                                               params,
                                               parallel=False,
                                               scale=1e6)
N_coalescence, R_coalescence = stabilizeCommunities(coalescence_plate)



### PLOT

# stack plot
fig, ax = myStackPlot(pd.concat([N_resident,N_invasive],axis=1))
ax.set_xlim(-0.5,1.5)
ax.set_xticklabels(['Resident\ncommunity','Invasive\ncommunity'])
ax.set_aspect(1.0/ax.get_data_ratio())

# dominant invasive fractions
f_singleinv = N_singleinv/N_singleinv.sum()
f_coalescence = N_coalescence/N_coalescence.sum()
f = [f_singleinv.loc[('F0','S2'),'W0'], f_coalescence.loc[('F0','S2'),'W0']]
fig, ax = plt.subplots()
ax.bar([0,1],f,width=0.5,color='black')
ax.set_xlim(-0.5,1.5)
ax.set_xticks([0,1])
ax.set_xticklabels(['Dominant\ninvading\nalone','Dominant\ninvading\nwith secondary'])
ax.set_aspect(1.0/ax.get_data_ratio())
ax.set_ylabel('Fraction of dominant invasive species\nin coalescence')
ax.set_ylim(0,0.5)
ax.set_yticks([0,0.2,0.4])
