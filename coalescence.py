# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:28:37 2020
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

'''
# set random and numpy's seed for reproducibility
myseed = 1000000 # 99
np.random.seed(myseed)
random.seed(myseed)
'''

import time
start_time = time.time()

# get community hierarchies? (this generates a complete data set but takes time)
getHierarchies = True



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

# rank plot
def rankPlot(N):
    
    # abundance threshold
    u = 1e-4
    
    # check number of wells
    n_wells = N.shape[1]
    
    # normalize to get fractions
    N = N/N.sum()
        
    # initialize figure
    fig, ax = plt.subplots(1)
    
    # plot
    for i in range(n_wells):
        
        # get ith column of N, remove elements below threshold and sort it
        Ni = N.iloc[:,i].to_numpy()
        Ni = Ni[Ni>u]
        Ni.sort()
        Ni = Ni[::-1]
        
        # plot
        ax.plot(Ni)
    
    ax.set_yscale("log")
    ax.set_ylabel('Relative abundance')
    ax.set_xlabel("Rank")
    ax.set_ylim(1e-4,1e0)
    ax.set_aspect(5)
        
    return([fig,ax])

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
    
    # position in the x axis
    xpos = np.arange(n_wells).tolist()
    xticklabs = np.add(xpos,[1]).tolist()
    
    # define cumulative sum of abundances
    cumsum = np.zeros(n_wells).tolist()
    
    # initialize figure
    fig, ax = plt.subplots(1)
    
    # plot bars - iterate across species and stack
    for i in range(n_species):
        ax.bar(xpos,N[i],color=cmap[i],bottom=cumsum)
        cumsum = np.add(N[i],cumsum).tolist()
        
    ax.set_xticks(xpos)
    ax.set_xticklabels(xticklabs)
    ax.set_ylabel('Fraction')
    ax.set_ylim(0,1)
        
    return [fig,ax]

# Bray-Curtis similarity
# p and q are lists (of absolute abundances, sum need not be 1)
def bray_curtis(p,q):

    # C_pq: sum of lesser values of common species
    c = sum([min(p[i],q[i]) for i in range(len(p)) if p[i]>0 and q[i]>0])
    
    # Bray-Curtis similarity
    bc = 2*c/(sum(p) + sum(q))
    return bc

# Kullback-Leibler divergence - using log2 so result is within [0,1]
# p and q are lists
def kldiv(p,q):
        
    # remove elements that are zero in both p and q (unused dimensions)
    n = [i for i in range(len(p)) if p[i]>0 or q[i]>0]
    p = [p[i] for i in n]
    q = [q[i] for i in n]
    
    # define this operation to handle the p[i]->0 limit
    def mylog(a,b):
        if a==0:
            return(0)
        else:
            return(a*math.log2(a/b))
    
    d = sum([mylog(p[i],q[i]) for i in range(len(p))])
    return d

# Jensen-Shannon similarity (1 - distance)
# p and q are lists
def jensen_shannon(p,q):
    
    # return an error if...
    if (not len(p)==len(q) or
    not math.isclose(sum(p),1) or
    not math.isclose(sum(q),1)):
        return "error"
    
    else:
        m = [(p[i]+q[i])/2 for i in range(len(p))]
        js_div = 1/2 * kldiv(p,m) + 1/2 * kldiv(q,m)
        js_dist = math.sqrt(js_div)
        return 1 - js_dist
    
# Jaccard similarity
# p and q are lists
def jaccard(p,q):
    
    intersection = sum([1 for i in range(len(p)) if p[i]>0 and q[i]>0])
    union = sum([1 for i in range(len(p)) if p[i]>0 or q[i]>0])
    
    return intersection/union

# endemic species overlap
def endemic(p,q,pq):
    
    overlap_p_pq = sum([1 for i in range(len(p))
                        if p[i]>0 and q[i]==0 and pq[i]>0])/sum([1 for i in range(len(p))
                                                                 if p[i]>0 and q[i]==0])
    overlap_q_pq = sum([1 for i in range(len(q))
                        if p[i]==0 and q[i]>0 and pq[i]>0])/sum([1 for i in range(len(q))
                                                                 if p[i]==0 and q[i]>0])
    
    return overlap_p_pq/(overlap_p_pq + overlap_q_pq)

# similarity index
# p: invasive
# q: resident
# pq: coalesced
def mysim(p,q,pq):
    
    # normalization
    p = [p_i/sum(p) for p_i in p]
    q = [q_i/sum(q) for q_i in q]
    pq = [pq_i/sum(pq) for pq_i in pq]
    
    sim = {'bray_curtis':bray_curtis(p,pq)/(bray_curtis(p,pq) + bray_curtis(q, pq)),
           'jensen_shannon':jensen_shannon(p,pq)/(jensen_shannon(p,pq) + jensen_shannon(q, pq)),
           'jaccard':jaccard(p,pq)/(jaccard(p,pq) + jaccard(q,pq)),
           'endemic':endemic(p,q,pq)}
                                                  
    return sim



# %%
### MODEL DEFINITION 

# general assumptions
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 100 # number of communities
assumptions['S'] = 50 # number of species sampled at initialization
assumptions['SA'] = [800, 800, 800] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 240 # 30 # number of generalists
assumptions['MA'] = [10, 10, 10] # [30, 30, 30] # number of resources per resource class
assumptions['l'] = 0.8 # leakage fraction

assumptions['response'] = 'type I'
assumptions['regulation'] = 'independent' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
assumptions['sampling'] = 'Gamma' # 'Binary', 'Gaussian', 'Uniform' or 'Gamma' (sampling of the matrix c)
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)

assumptions['q'] = 0.9 #0.9 # preference strength (0 for generalist and 1 for specialist)
assumptions['c0'] = 0.0 # background consumption rate in binary model
assumptions['c1'] = 1.0 # specific consumption rate in binary model
assumptions['sigc'] = 3 #3 # standard deviation of sum of consumption rates for Gaussian and Gamma models

assumptions['sparsity'] = 0.9 #0.05 # variability in secretion fluxes among resources (must be less than 1)  
assumptions['fs'] = 0.45 #0.45 # fraction of secretion flux to resources of the same type as the consumed one
assumptions['fw'] = 0.45 #0.45 # fraction of secretion flux to waste resources
assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['rs'] = 0.0 # control parameter for resource secretion: 0 means random secretions, 1 means species only secrete resources they can consume (relevant only if 'metabolism' is 'specific')

# parameters
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

'''
# plot initial state (10 communities max)
fig, ax = myStackPlot(resident_plate.N.iloc[:,0:min(10,resident_plate.N.shape[1])])
'''

# stabilize plates
N_resident, R_resident = stabilizeCommunities(resident_plate)
N_invasive, R_invasive = stabilizeCommunities(invasive_plate)

'''
# plot community composition after stabilization (10 communities max)
fig, ax = myStackPlot(resident_plate.N.iloc[:,0:min(10,resident_plate.N.shape[1])])
'''

'''
# plot ranks
rankPlot(N_resident)
'''

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
### DOMINANT SPECIES INVADING ALONE

# make plate (added dilution factor)
N0_singleinv = (resident_plate.N + monoculture_invasive_plate.N)/2*(1/100)
init_state_singleinv = (N0_singleinv,
                        init_state[1])
singleinv_plate = community_simulator.Community(init_state_singleinv,
                                                dynamics,
                                                params,
                                                parallel=par,
                                                scale=1e6)

# stabilize
N_singleinv, R_singleinv = stabilizeCommunities(singleinv_plate)



# %%
### HIERARCHY INDEX: HOW MUCH OF THE DOMINANTS' GROWTH COMES FROM THE PRIMARY RESOURCE?

# define dNdt on resource R0 only
def dNdt_R0(N,R,params):
    return community_simulator.usertools.MakeConsumerDynamics_R0(assumptions)(N,R,params)

# define function to calculate hierarchy (requires integration of dynamic equations; passed plate should be stabilized)
def communityHierarchy(plate,well):

    # dilute consumer and resource abundance of i-th well, replenish primary resource
    N = plate.N.values[:,well].copy() * (1/100)
    R = plate.R.values[:,well].copy() * (1/100)
    R[0] = R[0] + assumptions['R0_food']
    
    # keep track of species growth *on resource R0 only*
    N_R0 = N.copy()
    
    # find dominant
    dom = np.argmax(N)
    
    # integrate until stabilization (max_rel_change less than 1%)
    max_rel_change = np.inf
    dt = 1/1000 # time step for integration
    while max_rel_change > 0.1:
        
        # integrate species and resources abundance in the interval (t, t+dt)
        Nnew = N + dt*dNdt(N,R,plate.params)
        Nnew_R0 = N_R0 + dt*dNdt_R0(N,R,plate.params)
        Rnew = R + dt*dRdt(N,R,plate.params)
        
        # maximum relative change in abundance (omitting N=0)
        max_rel_change = np.max(np.divide(np.abs(Nnew - N),N,out=np.zeros_like(N),where=N!=0))
        
        # update species and resources abundances
        N = Nnew.copy()
        N_R0 = Nnew_R0.copy()
        R = Rnew.copy()
    
    # return hierarchy (fraction of dominant biomass that comes from resource R0)
    return N_R0[dom]/N[dom]

# calculate hierarchies of all communities (if applicable)
h_inv = ['NA' for i in range(assumptions['n_wells'])]
h_res = ['NA' for i in range(assumptions['n_wells'])]

if getHierarchies:
    print('Calculating community hierarchies...')
    for i in range(assumptions['n_wells']):
        print('... processing ' + '%s' % (i+1) + ' out of %s' % assumptions['n_wells'])
        h_inv[i] = communityHierarchy(invasive_plate,well=i)
        h_res[i] = communityHierarchy(resident_plate,well=i)



# %%
### STATISTICS

# community similarity metrics
s0 = ['NA' for i in range(assumptions['n_wells'])]
stats = {'bray_curtis':s0.copy(),
         'jensen_shannon':s0.copy(),
         'jaccard':s0.copy(),
         'endemic':s0.copy()}

for s in stats.keys():
    for i in range(assumptions['n_wells']):
        
        # community compositions as lists
        R = N_resident.iloc[:,i].tolist()
        I = N_invasive.iloc[:,i].tolist()
        C = N_coalescence.iloc[:,i].tolist()
        
        stats[s][i] = mysim(I,R,C)[s]
    
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
    
# get fraction in coalescence
f_coalescence = ['NA' for i in range(assumptions['n_wells'])]
f_coalescence_null = ['NA' for i in range(assumptions['n_wells'])]
F = N_coalescence/N_coalescence.sum()
F_null = N0_coalescence/N0_coalescence.sum()
F = F.droplevel(0)
F_null = F_null.droplevel(0)
for i in range(assumptions['n_wells']):
    f_coalescence[i] = F.loc[dominants_invasive[i]][i]
    f_coalescence_null[i] = F_null.loc[dominants_invasive[i]][i]
    
# get fracion when invading alone
f_singleinv = ['NA' for i in range(assumptions['n_wells'])]
f_singleinv_null = ['NA' for i in range(assumptions['n_wells'])]
F = N_singleinv/N_singleinv.sum()
F_null = N0_singleinv/N0_singleinv.sum()
F = F.droplevel(0)
F_null = F_null.droplevel(0)
for i in range(assumptions['n_wells']):
    f_singleinv[i] = F.loc[dominants_invasive[i]][i]
    f_singleinv_null[i] = F_null.loc[dominants_invasive[i]][i]



# %%
### SAVE DATA

# final data frame
data_out = pd.DataFrame(data={'f_pairwise':f_pairwise,
                              'q_bray_curtis':stats['bray_curtis'],
                              'q_jensen_shannon':stats['jensen_shannon'],
                              'q_jaccard':stats['jaccard'],
                              'q_endemic':stats['endemic'],
                              'f_singleinv':f_singleinv,
                              'f_coalescence':f_coalescence,
                              'h_inv':h_inv,
                              'h_res':h_res})
data_out.to_csv(os.path.join('.','data','simul_data.txt'),
                header=True,index=None,sep='\t',na_rep='NA')



# %%
### PLOTS

# Q vs fraction in pairwise competition
def q_vs_fraction():
    
    x = f_pairwise
    y = stats['bray_curtis']
    
    ### remove nan elements (this happens if pairwise competition ends up in extinction of both, need to investigate how this happens)
    ### SOLVED: this happens when none of the two competing species are able to grow over the dilution factor during the propagation time.
    ### It is an unusual case but it is possible and it does not imply code malfunction.
    n = [ni for ni in range(len(x)) if ~np.isnan(x[ni]) and ~np.isnan(y[ni])]
    x = [x[i] for i in n]
    y = [y[i] for i in n]
    
    fig, ax = plt.subplots()
    ax.scatter(x,y,edgecolors="black",facecolors="none")
    
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    dx = (max(x) - min(x))/10
    ax.plot([min(x)-dx,max(x)+dx],p([min(x)-dx,max(x)+dx]),c="black")
    
    ax.set_ylabel("Q\nCoalesced - Invasive")
    ax.set_xlabel("Frequency of invasive dominant species\nin pairwise competition")
    ax.set_xlim(-0.05,1.05)
    ax.set_ylim(-0.05,1.05)
    ax.set_aspect("equal")
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    
    fig

# fraction when invading alone vs with cohort
def cohort_vs_alone():
    
    x = f_singleinv
    y = f_coalescence
    
    fig, ax = plt.subplots()
    
    # filled regions
    w = 0.2 # width of the gray area (diagonal), must be between 0 and 1
    ax.fill([-0.1,0.1,0.1,-0.1],[-0.1,0.1,1.1,1.1],facecolor="green",alpha=0.15)
    ax.fill([-0.1,0.1,1.1,1.1],[-0.1,0.1,0.1,-0.1],facecolor="red",alpha=0.15)
    ax.fill([0,1.1*(1-w),1.1,1.1],[0,1.1,1.1,1.1*(1-w)],facecolor=[0.85,0.85,0.85])
    
    # scatterplot
    ax.scatter(x,y,edgecolors="black",facecolors="none",zorder=5)
    ax.plot([-0.1,1.1],[-0.1,1.1],'--k')
    
    ax.set_ylabel("Frequency of dominant species\ninvading with cohort")
    ax.set_xlabel("Frequency of dominant species\ninvading alone")
    ax.set_xlim(-0.05,1.05)
    ax.set_ylim(-0.05,1.05)
    ax.set_aspect("equal")
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    
    fig
    
# Q vs fraction in pairwise competition, split by strength of bottom-up co-selection
def q_vs_fraction_split():
    
    x = f_pairwise
    y = stats['bray_curtis']  
    
    ### remove nan elements (this happens if pairwise competition ends up in extinction of both, need to investigate how this happens)
    ### SOLVED: this happens when none of the two competing species are able to grow over the dilution factor during the propagation time.
    ### It is an unusual case but it is possible and it does not imply code malfunction.
    n = [ni for ni in range(len(x)) if ~np.isnan(x[ni]) and ~np.isnan(y[ni])]
    x = [x[i] for i in n]
    y = [y[i] for i in n]
    
    # split data
    x_buc = [x[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.05]
    y_buc = [y[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.05]
    
    x_nobuc = [x[i] for i in range(len(x)) if f_singleinv[i]>0.1 or f_coalescence[i]<0.05]
    y_nobuc = [y[i] for i in range(len(x)) if f_singleinv[i]>0.1 or f_coalescence[i]<0.05]
    
    # start plot
    fig, (ax1,ax2) = plt.subplots(1, 2)
    
    # positive bottom-up co-selection
    ax1.scatter(x_buc,y_buc,edgecolors="black",facecolors="none")
    z = np.polyfit(x_buc, y_buc, 1)
    p = np.poly1d(z)
    dx = (max(x_buc) - min(x_buc))/10
    ax1.plot([min(x_buc)-dx,max(x_buc)+dx],p([min(x_buc)-dx,max(x_buc)+dx]),c="black")
    
    ax1.set_ylabel("Q\nCoalesced - Invasive")
    #ax1.set_xlabel("Frequency of invasive dominant species\nin pairwise competition")
    ax1.set_xlim(-0.05,1.05)
    ax1.set_ylim(-0.05,1.05)
    ax1.set_aspect("equal")
    ax1.set_xticks([0,0.5,1])
    ax1.set_yticks([0,0.5,1])
    
    # no bottom-up co-selection
    ax2.scatter(x_nobuc,y_nobuc,edgecolors="black",facecolors="none")
    z = np.polyfit(x_nobuc, y_nobuc, 1)
    p = np.poly1d(z)
    dx = (max(x_nobuc) - min(x_nobuc))/10
    ax2.plot([min(x_nobuc)-dx,max(x_nobuc)+dx],p([min(x_nobuc)-dx,max(x_nobuc)+dx]),c="black")
    
    ax2.set_ylabel("")
    #ax2.set_xlabel("Frequency of invasive dominant species\nin pairwise competition")
    ax2.set_xlim(-0.05,1.05)
    ax2.set_ylim(-0.05,1.05)
    ax2.set_aspect("equal")
    ax2.set_xticks([0,0.5,1])
    ax2.set_yticks([0,0.5,1])
    
    # x axis label
    fig.text(0.5,0.04,"Frequency of invasive dominant species\nin pairwise competition",ha='center')
    
    # display plot
    fig

# make plots
q_vs_fraction()
cohort_vs_alone()
q_vs_fraction_split()



# %%
### TESTING AREA (keep commented when running from bash)

# execution time
print("%s s" % (time.time() - start_time))













