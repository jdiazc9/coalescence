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

# set random and numpy's seed manually for reproducibility
np.random.seed(0)
random.seed(0)

import time
start_time = time.time()



### USER-DEFINED FUNCTIONS

# generate a list of metabolic matrices
def MakeMatrices_multi(assumptions,mode='random'):
    if mode == 'random':
        return community_simulator.usertools.MakeMatrices(assumptions)[0], [community_simulator.usertools.MakeMatrices(assumptions)[1] for i in range(sum(assumptions['SA'])+assumptions['Sgen'])]
    elif mode == 'secrete_consumables':
        c, M = MakeMatrices_multi(assumptions,mode='random') # base case
        
        n = c == assumptions['c0']
        for i in range(len(M)): # loop through species
            ni = np.where(n.iloc[i,:])[0]
            M[i].iloc[ni,:] = 0 # M[i].iloc[ni,ni] = 0 also removes (sets to 0) the columns (in addition to the rows) corresponding to resources that the species i cannot consume, which should have no effect in practice
            M[i] = M[i]/M[i].sum()
        
        return c, M

# make resource dynamics if each species has its own metabolic matrix
def MakeResourceDynamics_multi(assumptions):
    """
    This function is adapted so that it returns a function (dRdt) which,
    in turn, ALWAYS expects params['D'] to be a LIST of metabolic matrices D_i
    of length equal to the total number of species (Stot).
    The i-th element of the list should the metabolic matrix of species i.
    
    Original description from the Community Simulator package:
    Construct resource dynamics. 'assumptions' must be a dictionary containing at least
    three entries:
    
    response = {'type I', 'type II', 'type III'} specifies nonlinearity of growth law
    
    regulation = {'independent','energy','mass'} allows microbes to adjust uptake
        rates to favor the most abundant accessible resources (measured either by
        energy or mass)
    
    supply = {'off','external','self-renewing'} sets choice of
        intrinsic resource dynamics
        
    Returns a function of N, R, and the model parameters, which itself returns the
        vector of resource rates of change dR/dt
    """
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
    
    """
    ### original in community-simulator package (single D matrix)
    J_out = lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)
    ###
    """
    
    """
    ### my first attempt at multiple matrices - very unefficient and slow bc it requires a lot of .dot() operations between large matrices
    def J_out(R,params):
        # get J_in so J_in(R,params) does not have to be re-evaluated for every element of the list
        Jin = J_in(R,params)
        # J_out is now a list of J_out_i matrices
        # J_out_i is the result of the operation (l*J_in)*D_i'
        # where D_i is the metabolic matrix for species i and ' indicates transposition
        #J_out = pd.Series(params['D']).apply(lambda x: (params['l']*J_in(R,params)).dot(x.T)).tolist()
        J_out = pd.Series(params['D']).apply(lambda x: (params['l']*Jin).dot(x.T)).tolist()
        # the final J_out comes from combining all the J_out_i row by row
        J_out = [J_out[i][i,:] for i in range(J_out[0].shape[0])]
        J_out = np.array(J_out)
        return J_out
    ###
    """
    
    ### my second attempt: using list comprehension and trimming J_in before .dot() to speed up
    J_out = lambda R,params: params['l']*np.array([ J_in(R,params)[i,:].dot(params['D'][i].T) for i in range(len(params['D'])) ])
    
    return lambda N,R,params: (h[assumptions['supply']](R,params)
                               -(J_in(R,params)/params['w']).T.dot(N)
                               +(J_out(R,params)/params['w']).T.dot(N))
    ###

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
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 10
assumptions['S'] = 100 # number of species sampled at initialization
assumptions['SA'] = [100, 100, 100] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 20 # number of generalists
assumptions['l'] = 0.8 # leakage fraction
assumptions['MA'] = [10, 10, 10] # [30, 30, 30] # number of resources per resource class
assumptions['response'] = 'type I'
assumptions['regulation'] = 'energy' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
assumptions['sampling'] = 'Binary' # 'Binary' or 'Gamma' (sampling of the matrix c)
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)
assumptions['c0'] = 0.0 # background consumption rate in binary model
assumptions['c1'] = 1 # specific consumption rate in binary model
assumptions['sparsity'] = 0.05 #0.05 # variability in secretion fluxes among resources (must be less than 1)  
assumptions['q'] = 0.9 #0.9 # preference strength (0 for generalist and 1 for specialist)
assumptions['metabolism'] = 'common' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['rs'] = 0 # resource secretion: 0 means random secretions, 1 means secretions are correlated to the species ability to consume the secreted resource (only used if 'metabolism' is 'specific')

# parameters
params = community_simulator.usertools.MakeParams(assumptions)

# if ''metabolism' is 'specific', make the list of D matrices...
if assumptions['metabolism'] == 'specific':
    params['c'], params['D'] = MakeMatrices_multi(assumptions,mode='random') # modes: 'random' generates random matrices according to assumptions, 'secrete_consumables' makes it so each species only secretes resources that be consumed by itself

# ...and make the dynamics depending on whether there is one or multiple D matrices
def dNdt(N,R,params):
    return community_simulator.usertools.MakeConsumerDynamics(assumptions)(N,R,params)

if assumptions['metabolism'] == 'specific':
    def dRdt(N,R,params):
        return MakeResourceDynamics_multi(assumptions)(N,R,params)
else:     
    def dRdt(N,R,params):
        return community_simulator.usertools.MakeResourceDynamics(assumptions)(N,R,params)

dynamics = [dNdt,dRdt]

# plot matrices c and D
fig,ax=plt.subplots()
sns.heatmap(params['c'],cmap='Greys',vmin=0,square=True,xticklabels=False,yticklabels=False,cbar=False,ax=ax)
ax.set_title('consumer matrix c')
fig
if isinstance(params['D'],list):
    Dplot = params['D'][4]
else:
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
                                               parallel=False,
                                               scale=1e6)

# initialize invasive plate
N0_invasive = init_state[0]
N0_invasive.loc[:,:] = sampleFromPool(1)
init_state_invasive = (N0_invasive,
                       init_state[1])
invasive_plate = community_simulator.Community(init_state_invasive,
                                               dynamics,
                                               params,
                                               parallel=False,
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
                                                  parallel=False,
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
                                                           parallel=False,
                                                           scale=1e6)

init_state_monoculture_invasive = (N0_monoculture_invasive,
                                   init_state_monoculture[1])
monoculture_invasive_plate = community_simulator.Community(init_state_monoculture_invasive,
                                                           dynamics,
                                                           params,
                                                           parallel=False,
                                                           scale=1e6)

# stabilize monocultures
monoculture_resident_plate.Propagate(1,compress_resources=False,compress_species=True)
monoculture_invasive_plate.Propagate(1,compress_resources=False,compress_species=True)



### PAIRWISE COMPETITION

# make plate (added dilution factor)
N0_pairwise = (monoculture_resident_plate.N + monoculture_invasive_plate.N)/2*(1/100)
init_state_pairwise = (N0_pairwise,
                       init_state[1])
pairwise_plate = community_simulator.Community(init_state_pairwise,
                                               dynamics,
                                               params,
                                               parallel=False,
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
                                                parallel=False,
                                                scale=1e6)

# stabilize
N_singleinv, R_singleinv = stabilizeCommunities(singleinv_plate)



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
    


### PLOTS

# Q vs fraction in pairwise competition
def q_vs_fraction():
    x = f_pairwise
    y = Q
    
    fig, ax = plt.subplots()
    ax.scatter(x,y,c="black")
    
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    dx = (max(x) - min(x))/10
    ax.plot([min(x)-dx,max(x)+dx],p([min(x)-dx,max(x)+dx]),c="black")
    
    ax.set_ylabel("Q\nCoalesced - Invasive")
    ax.set_xlabel("Frequency of invasive dominant species\nin pairwise competition")
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    
    fig

# Q vs fraction in pairwise competition
def cohort_vs_alone():
    x = f_singleinv
    y = f_coalescence
    
    fig, ax = plt.subplots()
    ax.scatter(x,y,c="black")
    ax.plot([0,1],[0,1],'--k')
    
    ax.set_ylabel("Frequency of dominant species\ninvading with cohort")
    ax.set_xlabel("Frequency of dominant species\ninvading alone")
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    
    fig

# make plots
q_vs_fraction()
cohort_vs_alone()



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

# test the behavior of the rs control parameter
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 10
assumptions['S'] = 100 # number of species sampled at initialization
assumptions['SA'] = [100, 100, 100] # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 20 # number of generalists
assumptions['l'] = 0.8 # leakage fraction
assumptions['MA'] = [10, 10, 10] # [30, 30, 30] # number of resources per resource class
assumptions['response'] = 'type I'
assumptions['regulation'] = 'energy' # 'independent' is standard, 'energy' or 'mass' enable preferential consumption of resources
assumptions['sampling'] = 'Binary' # 'Binary' or 'Gamma' (sampling of the matrix c)
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)
assumptions['c0'] = 0.0 # background consumption rate in binary model
assumptions['c1'] = 1 # specific consumption rate in binary model
assumptions['sparsity'] = 0.05 #0.05 # variability in secretion fluxes among resources (must be less than 1)  
assumptions['q'] = 0.9 #0.9 # preference strength (0 for generalist and 1 for specialist)
assumptions['metabolism'] = 'specific' # 'common' uses a common D matrix for all species, 'specific' uses a different matrix D for each species
assumptions['rs'] = 0 # resource secretion: 0 means random secretions, 1 means secretions are correlated to the species ability to consume the secreted resource (only used if 'metabolism' is 'specific')

c, D = community_simulator.usertools.myMakeMatrices(assumptions)
Dplot = D[0]
sns.heatmap(Dplot,cmap='Greys',vmin=0,square=True,xticklabels=False,yticklabels=False,cbar=False,ax=ax)
ax.set_title('metabolic matrix D')
fig















