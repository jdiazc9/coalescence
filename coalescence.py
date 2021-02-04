# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:28:37 2020

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
# import pandas as pd
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
    bc = 1 - 2*c/(sum(p) + sum(q))
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

# Jensen-Shannon distance
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
        return js_dist



### SEED AND STABILIZE COMMUNITIES

# assumptions

"""
# defaults for this simulations (these work)
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 20 # number of communities
assumptions['S'] = 100 # number of species sampled at initialization
assumptions['SA'] = [100, 100, 100] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 30 # number of generalists
assumptions['l'] = 0.8 # leakage fraction
assumptions['MA'] = [10, 10, 10] # [30, 30, 30] # number of resources per resource class

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

assumptions['sparsity'] = 0.6 #0.05 # variability in secretion fluxes among resources (must be less than 1)  
assumptions['fs'] = 0.45 #0.45 # fraction of secretion flux to resources of the same type as the consumed one
assumptions['fw'] = 0.45 #0.45 # fraction of secretion flux to waste resources
assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['rs'] = 0.0 # control parameter for resource secretion: 0 means random secretions, 1 means species only secrete resources they can consume (relevant only if 'metabolism' is 'specific')
"""

# try these assumptions
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 500 # number of communities
assumptions['S'] = 50 # number of species sampled at initialization
assumptions['SA'] = [200, 200, 200] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 60 # 30 # number of generalists
assumptions['l'] = 0.8 # leakage fraction
assumptions['MA'] = [10, 10, 10] # [30, 30, 30] # number of resources per resource class

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

assumptions['sparsity'] = 0.6 #0.05 # variability in secretion fluxes among resources (must be less than 1)  
assumptions['fs'] = 0.45 #0.45 # fraction of secretion flux to resources of the same type as the consumed one
assumptions['fw'] = 0.45 #0.45 # fraction of secretion flux to waste resources
assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['rs'] = 0.0 # control parameter for resource secretion: 0 means random secretions, 1 means species only secrete resources they can consume (relevant only if 'metabolism' is 'specific')

print(assumptions)

# parameters
params = community_simulator.usertools.MakeParams(assumptions)

# dynamics
def dNdt(N,R,params):
    return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)
def dRdt(N,R,params):
    return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)

dynamics = [dNdt,dRdt]

# plot matrices c and D
fig,ax=plt.subplots()
sns.heatmap(params['c'],cmap='Greys',vmin=0,square=True,xticklabels=False,yticklabels=False,cbar=False,ax=ax)
ax.set_title('consumer matrix c')
fig
fig,ax=plt.subplots()
sns.heatmap(params['c'].iloc[0:20,:],cmap='Greys',vmin=0,square=True,xticklabels=False,yticklabels=False,cbar=False,ax=ax)
ax.set_title('consumer matrix c (detail)')
fig
if assumptions['metabolism'] == 'specific':
    Dplot = params['D'][4]
elif assumptions['metabolism'] == 'common':
    Dplot = params['D']
fig,ax=plt.subplots()
sns.heatmap(Dplot,cmap='Greys',vmin=0,square=True,xticklabels=False,yticklabels=False,cbar=False,ax=ax)
ax.set_title('metabolic matrix D')
fig

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

# plot initial state
fig, ax = myStackPlot(resident_plate.N)

# stabilize plates
N_resident, R_resident = stabilizeCommunities(resident_plate)
N_invasive, R_invasive = stabilizeCommunities(invasive_plate)

# plot communities after stabilization
fig, ax = myStackPlot(N_resident)

# plot species abundance histogram
# fig, ax = plt.subplots()
# ax.hist(np.log10(N_resident+1e-20),bins=30,histtype='bar',stacked='True')
# ax.set_xlim(-12,4)
# ax.set_ylim(0,1000)
# fig

# plot ranks
rankPlot(N_resident)

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

# stabilize monocultures
monoculture_resident_plate.Propagate(1,compress_resources=False,compress_species=True)
monoculture_invasive_plate.Propagate(1,compress_resources=False,compress_species=True)



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



### COHORT INVADING WITHOUT DOMINANT SPECIES

# make cohort plate
invasive_cohorts_plate = invasive_plate.copy()
invasive_cohorts_plate.N = invasive_cohorts_plate.N.droplevel(0)
for i in range(assumptions['n_wells']):
    invasive_cohorts_plate.N.loc[dominants_invasive[i]][i] = 0
invasive_cohorts_plate.N.index = invasive_plate.N.index
    
# make plate (added dilution factor)
N0_cohortinv = (resident_plate.N + invasive_cohorts_plate.N)/2*(1/100)
init_state_cohortinv = (N0_cohortinv,
                        init_state[1])
cohortinv_plate = community_simulator.Community(init_state_cohortinv,
                                                dynamics,
                                                params,
                                                parallel=par,
                                                scale=1e6)

# stabilize
N_cohortinv, R_cohortinv = stabilizeCommunities(cohortinv_plate)



### GET STATISTICS

# similarity index (based on jensen_shannon, consider alternatives)
# note: jensen-shannon is a divergence, 
# sqrt(JS) is a distance
# 1-sqrt(JS) is a similarity measure
def mysim(p,q):
    p = [p_i/sum(p) for p_i in p]
    q = [q_i/sum(q) for q_i in q]
    return 1 - math.sqrt(jensen_shannon(p,q))

# get Q for each community
Q = ['NA' for i in range(assumptions['n_wells'])]
for i in range(assumptions['n_wells']):
    
    # community compositions as lists
    R = N_resident.iloc[:,i].tolist()
    I = N_invasive.iloc[:,i].tolist()
    C = N_coalescence.iloc[:,i].tolist()
    
    # i-th element of Q
    Q [i] = mysim(I,C)/(mysim(I,C) + mysim(R,C))
    
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
    
# get the ratio of the abundances of the dominant species when growing alone vs when growing with cohort
kratio_resident = ['NA' for i in range(assumptions['n_wells'])]
kratio_invasive = ['NA' for i in range(assumptions['n_wells'])]
N_resident = N_resident.droplevel(0)
N_invasive = N_invasive.droplevel(0)
N_resident_monoculture = monoculture_resident_plate.N.droplevel(0)
N_invasive_monoculture = monoculture_invasive_plate.N.droplevel(0)
for i in range(assumptions['n_wells']):
    kratio_resident[i] = N_resident.loc[dominants_resident[i]][i]/N_resident_monoculture.loc[dominants_resident[i]][i]
    kratio_invasive[i] = N_invasive.loc[dominants_invasive[i]][i]/N_invasive_monoculture.loc[dominants_invasive[i]][i]
kratio_resident = [math.log10(x) for x in kratio_resident]
kratio_invasive = [math.log10(x) for x in kratio_invasive]
kratio_diff = [kratio_invasive[i] - kratio_resident[i] for i in range(assumptions['n_wells'])]

# get cohort invasiveness (ci)
ci = ['NA' for i in range(assumptions['n_wells'])]
N_cohortinv = N_cohortinv.droplevel(0)
for i in range(assumptions['n_wells']):
    
    # identify species in invasive cohort
    cohort = N_invasive.index[N_invasive.iloc[:,i] > 0]
    cohort = [x for x in cohort if not(x==dominants_invasive[i])]
    
    # how many of those are in the final community? (resident + cohort invading alone)
    cohort_afterinv = N_cohortinv.iloc[:,i]
    cohort_afterinv = cohort_afterinv[cohort]
    
    # cohort invasivity is the fraction of the species in the cohort that were able to invade
    ci[i] = sum(cohort_afterinv > 0)/len(cohort)
    
# get communities bottom-up cohesiveness (buc)
buc_resident = ['NA' for i in range(assumptions['n_wells'])]
buc_invasive = ['NA' for i in range(assumptions['n_wells'])]

sigma = {'type I': lambda R,params: params['c']*R,
         'type II': lambda R,params: params['c']*R/(1+params['c']*R/params['sigma_max']),
         'type III': lambda R,params: (params['c']*R)**params['n']/(1+(params['c']*R)**params['n']/params['sigma_max'])
    }

u = {'independent': lambda x,params: 1.,
     'energy': lambda x,params: (((params['w']*x)**params['nreg']).T
                                  /np.sum((params['w']*x)**params['nreg'],axis=1)).T,
     'mass': lambda x,params: ((x**params['nreg']).T/np.sum(x**params['nreg'],axis=1)).T
    }

h = {'off': lambda R,params: 0.,
     'external': lambda R,params: (params['R0']-R)/params['tau'],
     'self-renewing': lambda R,params: params['r']*R*(params['R0']-R),
     'predator': lambda R,params: params['r']*R*(params['R0']-R)-params['u']*R
    }

J_in = lambda R,params: (u[assumptions['regulation']](params['c']*R,params)
                         *params['w']*sigma[assumptions['response']](R,params))
J_out = {'common': lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T),
         'specific': lambda R,params: params['l']*np.array([ J_in(R,params)[i,:].dot(params['D'][i].T) for i in range(len(params['D'])) ])
        }

# resident plate
plate = resident_plate.copy()
dominants = dominants_resident.copy()
for i in range(assumptions['n_wells']):
    
    # remove dominant
    plate.N.iloc[int(dominants[i].split('S',1)[1]),i] = 0
    
    # get secretions of dominant when feeding on the primary resource
    if assumptions['metabolism'] == 'specific':
        R0 = plate.params['D'][int(dominants[i].split('S',1)[1])]
        R0 = R0[:,0]
    elif assumptions['metabolism'] == 'common':
        R0 = plate.params['D'][:,0]
        
    # get secretions of cohort when feeding on the dominant's secretions
    resource_flux = (J_out[assumptions['metabolism']](R0,plate.params)/plate.params['w']).T.dot(plate.N.iloc[:,i])
    resource_flux = resource_flux.tolist()
    
    # get uptake rates of dominant species
    dominant_uptake_rates = plate.params['c'][int(dominants[i].split('S',1)[1]),:]
    dominant_uptake_rates = dominant_uptake_rates.tolist()
    
    # calculate overlap: use Jensen-Shannon
    ### FIXME: try other ways to quantify this overlap?
    x = [dominant_uptake_rates[j]/sum(dominant_uptake_rates) for j in range(len(dominant_uptake_rates))]
    y = [resource_flux[j]/sum(resource_flux) for j in range(len(resource_flux))]
    buc_resident[i] = 1 - math.sqrt(jensen_shannon(x,y))
    
# invasive plate
plate = invasive_plate.copy()
dominants = dominants_invasive.copy()
for i in range(assumptions['n_wells']):
    
    # remove dominant
    plate.N.iloc[int(dominants[i].split('S',1)[1]),i] = 0
    
    # get secretions of dominant when feeding on the primary resource
    if assumptions['metabolism'] == 'specific':
        R0 = plate.params['D'][int(dominants[i].split('S',1)[1])]
        R0 = R0[:,0]
    elif assumptions['metabolism'] == 'common':
        R0 = plate.params['D'][:,0]
        
    # get secretions of cohort when feeding on the dominant's secretions
    resource_flux = (J_out[assumptions['metabolism']](R0,plate.params)/plate.params['w']).T.dot(plate.N.iloc[:,i])
    resource_flux = resource_flux.tolist()
    
    # get uptake rates of dominant species
    dominant_uptake_rates = plate.params['c'][int(dominants[i].split('S',1)[1]),:]
    dominant_uptake_rates = dominant_uptake_rates.tolist()
    
    # calculate overlap: use Jensen-Shannon
    ### FIXME: try other ways to quantify this overlap?
    x = [dominant_uptake_rates[j]/sum(dominant_uptake_rates) for j in range(len(dominant_uptake_rates))]
    y = [resource_flux[j]/sum(resource_flux) for j in range(len(resource_flux))]
    buc_invasive[i] = 1 - math.sqrt(jensen_shannon(x,y))
   
# difference in bottom-up cohesiveness
buc_diff = [buc_invasive[j]-buc_resident[j] for j in range(assumptions['n_wells'])]
    


### PLOTS

# Q vs fraction in pairwise competition
def q_vs_fraction():
    
    x = f_pairwise
    y = Q
    
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
    
# k-ratio histogram
def kratio_hist():
    
    color_resident = [55/255,126/255,184/255]
    color_invasive = [255/255,127/255,0/255]
    
    x = np.linspace(min(kratio_resident+kratio_invasive),
                    max(kratio_resident+kratio_invasive),
                    num=50).tolist()
    
    fig, ax = plt.subplots()
    
    ax.hist(kratio_resident,
            bins=x,
            histtype='stepfilled',
            color=color_resident+[0.5],
            edgecolor=color_resident+[1])
    
    ax.hist(kratio_invasive,
            bins=x,
            histtype='stepfilled',
            color=color_invasive+[0.5],
            edgecolor=color_invasive+[1])
    
    ax.set_xlabel("log$_{10}$ K-ratio")
    ax.set_ylabel("Fraction")
    ax.set_aspect(1.0/ax.get_data_ratio()) # square axes even if different axes limits
    plt.legend(['Resident','Invasive'])
    
    fig
    
# difference in k-ratio vs. cohort invasiveness
def kratio_vs_ci():
    
    x = ci
    y = kratio_diff
    
    # identify dots in green area
    xg = [x[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.1]
    yg = [y[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.1]
    
    # identify dots outside of green area
    xo = [x[i] for i in range(len(x)) if not(f_singleinv[i]<0.1) or not(f_coalescence[i]>0.1)]
    yo = [y[i] for i in range(len(x)) if not(f_singleinv[i]<0.1) or not(f_coalescence[i]>0.1)]
    
    fig, ax = plt.subplots()
    
    ax.scatter(xo,yo,edgecolors=[0.5, 0.5, 0.5, 1],facecolors=[0.5, 0.5, 0.5, 0.5],zorder=1)
    ax.scatter(xg,yg,edgecolors=[0, 0.75, 0, 1],facecolors=[0, 0.75, 0, 0.5],zorder=1)
    
    ax.set_xlabel("Cohort invasiveness")
    ax.set_ylabel("log$_{10}$ K-ratio\nInvasive - Resident")
    ax.set_aspect(1.0/ax.get_data_ratio()) # square axes even if different axes limits
    
    fig
    
# difference in bottom-up cohesiveness vs. cohort invasiveness
def buc_vs_ci():
    
    x = ci
    y = buc_diff
    
    # identify dots in green area
    xg = [x[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.1]
    yg = [y[i] for i in range(len(x)) if f_singleinv[i]<0.1 and f_coalescence[i]>0.1]
    
    # identify dots outside of green area
    xo = [x[i] for i in range(len(x)) if not(f_singleinv[i]<0.1) or not(f_coalescence[i]>0.1)]
    yo = [y[i] for i in range(len(x)) if not(f_singleinv[i]<0.1) or not(f_coalescence[i]>0.1)]
    
    fig, ax = plt.subplots()
    
    ax.scatter(xo,yo,edgecolors=[0.5, 0.5, 0.5, 1],facecolors=[0.5, 0.5, 0.5, 0.5],zorder=1)
    ax.scatter(xg,yg,edgecolors=[0, 0.75, 0, 1],facecolors=[0, 0.75, 0, 0.5],zorder=1)
    
    ax.set_xlabel("Cohort invasiveness")
    ax.set_ylabel("Bottom-up cohesiveness\nInvasive - Resident")
    ax.set_aspect(1.0/ax.get_data_ratio()) # square axes even if different axes limits
    
    fig
    

# make plots
q_vs_fraction()
cohort_vs_alone()
kratio_hist()
kratio_vs_ci()
buc_vs_ci()



### TESTING AREA (keep commented when running from bash)

"""
# test how much time it takes to propagate a plate for T=1
import time
plate = resident_plate.copy()

start_time = time.time()
plate.Propagate(1)
print("%s s" % (time.time() - start_time))
"""

print("%s s" % (time.time() - start_time))













