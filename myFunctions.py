# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 15:55:43 2020

@author: Juan
"""

# libraries
import math

# Kullback-Leibler divergence - using log2 so result is within [0,1]
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
def js(p,q):
    
    # return an error if...
    if (not len(p)==len(q) or
    not sum(p)==1 or
    not sum(q)==1 or
    all(p_i==0 for p_i in p) or
    all(q_i==0 for q_i in q)):
        return "error"
    
    else:
        m = [(p[i]+q[i])/2 for i in range(len(p))]
        js_div = 1/2 * kldiv(p,m) + 1/2 * kldiv(q,m)
        js_dist = math.sqrt(js_div)
        return js_dist
