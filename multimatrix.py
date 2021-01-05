# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:53:14 2020

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
#import matplotlib.pyplot as plt
#import seaborn as sns
#colors = sns.color_palette()
# %matplotlib inline

import random
#import math
#import copy

# set random and numpy's seed manually for reproducibility
np.random.seed(0)
random.seed(0)





# assumptions
assumptions = community_simulator.usertools.a_default.copy()
assumptions['n_wells'] = 10
assumptions['S'] = 100 # number of species sampled at initialization
assumptions['SA'] = 200 # [100, 100, 100] # number of species per specialist family
assumptions['Sgen'] = 0 # number of generalists
assumptions['l'] = 0.8 # leakage fraction
assumptions['MA'] = 100 # [30, 30, 30] # number of resources per resource class
assumptions['response'] = 'type I'
assumptions['sampling'] = 'Binary' # 'Gamma'
assumptions['supply'] = 'off'
assumptions['R0_food'] = 1000
assumptions['m'] = 0 # turn off mortality (?)
assumptions['c0'] = 0.01 # background consumption rate in binary model
assumptions['c1'] = 1 # specific consumption rate in binary model

params = community_simulator.usertools.MakeParams(assumptions)
Dlist = [params['D'], params['D'], params['D']]

N, R = community_simulator.usertools.MakeInitialState(assumptions)
N = N.iloc[:,0]
R = R.iloc[:,0]





# taken from community_simulator.usertools as is
def MakeConsumerDynamics(assumptions):
    """
    Construct resource dynamics. 'assumptions' must be a dictionary containing at least
    three entries:
    
    response = {'type I', 'type II', 'type III'} specifies nonlinearity of growth law
    
    regulation = {'independent','energy','mass'} allows microbes to adjust uptake
        rates to favor the most abundant accessible resources (measured either by
        energy or mass)
    
    supply = {'off','external','self-renewing','predator'} sets choice of
        intrinsic resource dynamics
        
    Returns a function of N, R, and the model parameters, which itself returns the
        vector of consumer rates of change dN/dt
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
    
    J_in = lambda R,params: (u[assumptions['regulation']](params['c']*R,params)
                             *params['w']*sigma[assumptions['response']](R,params))
    J_growth = lambda R,params: (1-params['l'])*J_in(R,params)
    
    return lambda N,R,params: params['g']*N*(np.sum(J_growth(R,params),axis=1)-params['m'])

def MakeResourceDynamics(assumptions):
    """
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
    J_out = lambda R,params: (params['l']*J_in(R,params)).dot(params['D'].T)
    
    return lambda N,R,params: (h[assumptions['supply']](R,params)
                               -(J_in(R,params)/params['w']).T.dot(N)
                               +(J_out(R,params)/params['w']).T.dot(N))



###
# this function expects D to be a list of metabolic matrices (one per species)
def MakeResourceDynamics_multiMatrix(assumptions):
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
    def J_out(R,params):
        # get J_in so J_in(R,params) does not have to be re-evaluated for every element of the list
        Jin = J_in(R,params)
        # J_out is now a list of J_out_i matrices
        # J_out_i is the result of the operation (l*J_in)*D_i'
        # where D_i is the metabolic matrix for species i and ' indicates transposition
        #J_out = pd.Series(params['D']).apply(lambda x: (params['l']*J_in(R,params)).dot(x.T)).tolist()
        J_out = pd.Series(params['D']).apply(lambda x: (params['l']*Jin).dot(x.T)).tolist()
        # the final J_out comes from combining all the J_out_i row by row
        J_out = [J_out[i].iloc[i,:] for i in range(J_out[0].shape[0])]
        J_out = pd.concat(J_out,axis=1).T
        return J_out
    
    return lambda N,R,params: (h[assumptions['supply']](R,params)
                               -(J_in(R,params)/params['w']).T.dot(N)
                               +(J_out(R,params)/params['w']).T.dot(N))

