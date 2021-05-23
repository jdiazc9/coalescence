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
# import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette()
# %matplotlib inline

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
    dilution_factor = 0.01
    time_per_cycle = 1
    
    # get starting abundance of species in plate
    Nold = plate.N.copy()
    
    # run stabilization protocol for one transfer
    Ntraj, Rtraj = plate.RunExperiment(f = np.diag(np.ones(plate.n_wells))*dilution_factor,
                                       T = time_per_cycle,
                                       npass = 1,
                                       compress_resources=False,
                                       compress_species=True)
    
    # get new abundances
    Nnew = plate.N.copy()
    
    # maximum relative change in abundance
    d = np.max(np.abs(Nnew-Nold)/Nold)
    
    # run stabilization cycles until maximum relative change in abundance is below 1% (generational equlibrium achieved)
    while d[0] > 0.01:
        Nold = plate.N.copy()
        Ntraj, Rtraj = plate.RunExperiment(f = np.diag(np.ones(plate.n_wells))*dilution_factor,
                                       T = time_per_cycle,
                                       npass = 1,
                                       compress_resources=False,
                                       compress_species=True)
        Nnew = plate.N.copy()
        d = np.max(np.abs(Nnew-Nold)/Nold)
    
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

assumptions['sampling'] = 'Gamma'
assumptions['response'] = 'type I'
assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
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



# %%
### COALESCENCE EXPERIMENTS SETUP

# define experiment function (so it is iterable)
def competeMinimalCommunities(c_11 = 1.5, c_31 = 1, c_13 = 0, c_35 = 0, c_22 = 100, c_44 = 100):

    # consumer matrix
    params['c'].iloc[:,:] = np.array([[c_11,    0, c_13,    0,    0],
                                      [   0, c_22,    0,    0,    0],
                                      [c_31,    0,    0,    0, c_35],
                                      [   0,    0,    0, c_44,    0]])
    
    # metabolic matrices
    for i in range(len(params['D'])):
        params['D'][i].iloc[:,:] = np.identity(int(assumptions['MA']))
    
    params['D'][0].iloc[:,0] = np.array([0, 1, 0, 0, 0])
    params['D'][1].iloc[:,1] = np.array([0, 0, 1, 0, 0])
    params['D'][2].iloc[:,0] = np.array([0, 0, 0, 1, 0])
    params['D'][3].iloc[:,3] = np.array([0, 0, 0, 0, 1])
    
    # generic initial state
    init_state = community_simulator.usertools.MakeInitialState(assumptions)
    
    # resident community
    N0_res = init_state[0].copy()
    N0_res.iloc[:,0] = 0
    N0_res.iloc[0:2,0] = 1
    plate_res = community_simulator.Community((N0_res,init_state[1]),
                                              dynamics,
                                              params,
                                              parallel=par,
                                              scale=1e6)
    
    # invasive community
    N0_inv = init_state[0].copy()
    N0_inv.iloc[:,0] = 0
    N0_inv.iloc[2:5,0] = 1
    plate_inv = community_simulator.Community((N0_inv,init_state[1]),
                                              dynamics,
                                              params,
                                              parallel=par,
                                              scale=1e6)
    
    # stabilize primary communities
    N_res, R_res = stabilizeCommunities(plate_res)
    N_inv, R_inv = stabilizeCommunities(plate_inv)
    
    # coalesced community
    plate_coa = community_simulator.Community(((N_inv+N_res)/2,init_state[1]),
                                              dynamics,
                                              params,
                                              parallel=par,
                                              scale=1e6)
    
    # stabilize coalesced community
    N_coa, R_coa = stabilizeCommunities(plate_coa)
    
    # monoculture: resident dominant
    N0_res_mono = init_state[0].copy()
    N0_res_mono.iloc[:,0] = 0
    N0_res_mono.iloc[0,0] = 1
    plate_res_mono = community_simulator.Community((N0_res_mono,init_state[1]),
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # monoculture: invasive dominant
    N0_inv_mono = init_state[0].copy()
    N0_inv_mono.iloc[:,0] = 0
    N0_inv_mono.iloc[2,0] = 1
    plate_inv_mono = community_simulator.Community((N0_inv_mono,init_state[1]),
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # stabilize monocultures
    N_res_mono, R_res_mono = stabilizeCommunities(plate_res_mono)
    N_inv_mono, R_inv_mono = stabilizeCommunities(plate_inv_mono)
    
    # dominant invading alone
    N0_singleinv = N_res.copy()
    N0_singleinv.iloc[2,0] = N_inv_mono.iloc[2,0]
    plate_singleinv = community_simulator.Community((N0_singleinv,init_state[1]),
                                                    dynamics,
                                                    params,
                                                    parallel=par,
                                                    scale=1e6)
    
    # stabilize plate of dominant invading alone
    N_singleinv, R_singleinv = stabilizeCommunities(plate_singleinv)
    
    # pairwise competition of dominants
    N0_pairwise = N_res_mono + N_inv_mono
    plate_pairwise = community_simulator.Community((N0_pairwise,init_state[1]),
                                                   dynamics,
                                                   params,
                                                   parallel=par,
                                                   scale=1e6)
    
    # stabilize plate of pairwise competition
    N_pairwise, R_pairwise = stabilizeCommunities(plate_pairwise)
    
    # get statistics
    Q = mysim(N_inv.iloc[:,0].tolist(),N_res.iloc[:,0].tolist(),N_coa.iloc[:,0].tolist())
    f_pairwise = N_pairwise.iloc[2,0]/N_pairwise.iloc[:,0].sum()
    f_coalescence = N_coa.iloc[2,0]/N_coa.iloc[:,0].sum()
    f_singleinv = N_singleinv.iloc[2,0]/N_singleinv.iloc[:,0].sum()
    
    # return statistics
    return {'Q':Q,
            'f_pairwise':f_pairwise,
            'f_coalescence':f_coalescence,
            'f_singleinv':f_singleinv}



# %%
### PARAMETER GRIDS

# set resolution
n = 20 # number of breaks per axis (total experiments to run: n^2)


# top-down vs. top-down
c_11 = np.ones((n,n))
c_31 = np.tile(10**np.linspace(-0.5,0.5,num=n),(n,1))
c_13 = np.zeros((n,n))
c_35 = np.zeros((n,n))


'''
# top-down vs. bidirectional
c_11 = np.ones((n,n))
c_31 = np.tile(10**np.linspace(-0.5,0.5,num=n),(n,1))
c_13 = np.zeros((n,n))
c_35 = np.tile(10**np.linspace(3,-3,num=n).reshape((n,1)),(1,n))
'''

'''
# bidirectional vs. top-down
c_11 = np.ones((n,n))
c_31 = np.tile(10**np.linspace(-0.5,0.5,num=n),(n,1))
c_13 = np.tile(10**np.linspace(3,-3,num=n).reshape((n,1)),(1,n))
c_35 = np.zeros((n,n))
'''

'''
# bidirectional vs. bidirectional
c_11 = np.ones((n,n))
c_31 = np.tile(10**np.linspace(-0.5,0.5,num=n),(n,1))
c_13 = np.ones((n,n))
c_35 = np.tile(10**np.linspace(3,-3,num=n).reshape((n,1)),(1,n))
'''


# %%
### SWEEP PARAMETER SPACE

# initialize variables to store data
Q = [['NA' for i in range(n)] for i in range(n)]
f_pairwise = [['NA' for i in range(n)] for i in range(n)]
f_coalescence = [['NA' for i in range(n)] for i in range(n)]
f_singleinv = [['NA' for i in range(n)] for i in range(n)]

# run experiments iteratively through grid
for i in range(n):
    for j in range(n):
        exp_ij = competeMinimalCommunities(c_11 = c_11[i,j],
                                           c_31 = c_31[i,j],
                                           c_13 = c_13[i,j],
                                           c_35 = c_35[i,j])
        Q[i][j] = exp_ij['Q']
        f_pairwise[i][j] = exp_ij['f_pairwise']
        f_coalescence[i][j] = exp_ij['f_coalescence']
        f_singleinv[i][j] = exp_ij['f_singleinv']
        
        # check progress
        print("p_1: %s" % (i+1) + " out of " + str(n) + ", p_2: %s" % (j+1) + " out of " + str(n))
        print("%s seconds" % (time.time() - start_time))
        print('')



# %%
### PLOT Q MAP

from matplotlib.colors import LinearSegmentedColormap

# make colormap
myblue = np.array([94, 142, 194])/256
myyellow = np.array([232, 174, 34])/256
mycmap = LinearSegmentedColormap.from_list('mycmap',[myyellow,myblue],N=100)

# make plot
plt.pcolormesh(Q,edgecolors='none',cmap=mycmap)

# adjust axis
ax = plt.gca()
ax.set_xticks(np.array(range(n)) + 0.5)
ax.set_yticks(np.array(range(n)) + 0.5)
ax.set_xticklabels(np.array(range(n)))
ax.set_yticklabels(np.array(range(n)))
ax.set_aspect(1)
ax.invert_yaxis()

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.show()



# %%
### PLOT F(ALONE) MAP

# make plot
plt.pcolormesh(f_singleinv,edgecolors='none',cmap='binary')
plt.clim(-0.1,1.1)

# adjust axis
ax = plt.gca()
ax.set_xticks(np.array(range(n)) + 0.5)
ax.set_yticks(np.array(range(n)) + 0.5)
ax.set_xticklabels(np.array(range(n)))
ax.set_yticklabels(np.array(range(n)))
ax.set_aspect(1)
ax.invert_yaxis()

ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)



# %%
### PLOT F(PAIRWISE) VS. Q

x = np.array(f_pairwise).reshape(n**2)
y = np.array(Q).reshape(n**2)

fig, ax = plt.subplots()
ax.scatter(x,y,edgecolors="black",facecolors="none")

ax.set_ylabel("Q")
ax.set_xlabel("f (pairwise)")
ax.set_xlim(-0.05,1.05)
ax.set_ylim(-0.05,1.05)
ax.set_aspect("equal")
ax.set_xticks([0,0.5,1])
ax.set_yticks([0,0.5,1])

fig


# %%
### PLOT F(WITH COHORT) VS. F(ALONE)

x = np.array(f_singleinv).reshape(n**2)
y = np.array(f_coalescence).reshape(n**2)

fig, ax = plt.subplots()
ax.scatter(x,y,edgecolors="black",facecolors="none")

ax.set_ylabel("f (with cohort)")
ax.set_xlabel("f (alone)")
ax.set_xlim(-0.05,1.05)
ax.set_ylim(-0.05,1.05)
ax.set_aspect("equal")
ax.set_xticks([0,0.5,1])
ax.set_yticks([0,0.5,1])

fig


