# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:09:01 2020

@author: Juan
"""
# original versions of community-simulator methods
def CompressParams(not_extinct_consumers,not_extinct_resources,params,dimensions,S,M):
    params_comp = params.copy()
    if 'SxM' in dimensions.keys():
        for item in dimensions['SxM']:
            if item in params_comp.keys():
                assert np.shape(params_comp[item])==(S,M), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                params_comp[item]=params_comp[item][not_extinct_consumers,:]
                params_comp[item]=params_comp[item][:,not_extinct_resources]
    if 'MxM' in dimensions.keys():
        for item in dimensions['MxM']:
            if item in params_comp.keys():
                assert np.shape(params_comp[item])==(M,M), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                params_comp[item]=params_comp[item][not_extinct_resources,:]
                params_comp[item]=params_comp[item][:,not_extinct_resources]
    if 'SxS' in dimensions.keys():
        for item in dimensions['SxS']:
            if item in params_comp.keys():
                assert np.shape(params_comp[item])==(S,S), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                params_comp[item]=params_comp[item][not_extinct_consumers,:]
                params_comp[item]=params_comp[item][:, not_extinct_consumers]
    if 'S' in dimensions.keys():
        for item in dimensions['S']:
            if item in params_comp.keys():
                if type(params_comp[item]) == np.ndarray:
                    assert len(params_comp[item])==S, 'Invalid length for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                    params_comp[item]=params_comp[item][not_extinct_consumers]
    if 'M' in dimensions.keys():
        for item in dimensions['M']:
            if item in params_comp.keys():
                if type(params_comp[item]) == np.ndarray:
                    assert len(params_comp[item])==M, 'Invalid length for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                    params_comp[item]=params_comp[item][not_extinct_resources]  
    return params_comp


# original versions of community-simulator methods
def CompressParams(not_extinct_consumers,not_extinct_resources,params,dimensions,S,M):
    params_comp = params.copy()
    if 'SxM' in dimensions.keys():
        for item in dimensions['SxM']:
            if item in params_comp.keys():
                assert np.shape(params_comp[item])==(S,M), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                params_comp[item]=params_comp[item][not_extinct_consumers,:]
                params_comp[item]=params_comp[item][:,not_extinct_resources]
                
    if 'MxM' in dimensions.keys():
        for item in dimensions['MxM']:
            if item in params_comp.keys():
                ### FIXME: added this condition: is D a np.ndarray or a list?
                if isinstance(params_comp[item],np.ndarray): ### FIXME: if D is a pd.DataFrame, proceed noemally (compress it)
                    assert np.shape(params_comp[item])==(M,M), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                    params_comp[item]=params_comp[item][not_extinct_resources,:]
                    params_comp[item]=params_comp[item][:,not_extinct_resources]
                elif isinstance(params_comp[item],list): ### FIXME: if D is a list, compress every element of the list (keep the list structure to wrap them)
                    params_comp[item]=pd.Series(params_comp[item]).apply(lambda x: x[not_extinct_resources,:]).tolist()
                    params_comp[item]=pd.Series(params_comp[item]).apply(lambda x: x[:,not_extinct_resources]).tolist()
   
    if 'SxS' in dimensions.keys():
        for item in dimensions['SxS']:
            if item in params_comp.keys():
                assert np.shape(params_comp[item])==(S,S), 'Invalid shape for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                params_comp[item]=params_comp[item][not_extinct_consumers,:]
                params_comp[item]=params_comp[item][:, not_extinct_consumers]
    if 'S' in dimensions.keys():
        for item in dimensions['S']:
            if item in params_comp.keys():
                if type(params_comp[item]) == np.ndarray:
                    assert len(params_comp[item])==S, 'Invalid length for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                    params_comp[item]=params_comp[item][not_extinct_consumers]
    if 'M' in dimensions.keys():
        for item in dimensions['M']:
            if item in params_comp.keys():
                if type(params_comp[item]) == np.ndarray:
                    assert len(params_comp[item])==M, 'Invalid length for ' + item + '. Please update dimensions dictionary with correct dimensions.'
                    params_comp[item]=params_comp[item][not_extinct_resources]  
    return params_comp