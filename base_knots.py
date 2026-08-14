import numpy as np
import random
import pandas as pd
import gc

# IMPORTING A SUBSET OF BASE KNOTS

def import_base_knots(x):
    p05 = pd.read_csv(x)
    data_g = np.array(p05)

    return data_g

def generate_base_knots_dict(x, n=0, N=250, sf=20):
    '''generate a dictionary containing the specified range of knots
    
    Parameters:
    x: csv file storing the knot pts
    n, N: indeces to specify the range of knots to import
    sf: scaling factor multiplied by the coordddinates to avoid float point innacuracies during rasterisation

    Returns: dictionary containing 'knot id: coords array' pairs 
    '''
    # max N is 12965 in the knotinfo dataset
    f = import_base_knots(x)[n:N]
    pts_list=[]    
    knotid=[]
    base_knots1 = {}
    
    for k in range(f.shape[0]):      # coords
        a = f[k,3].split(';')
        b = []
        for n in range(len(a)):
            b.append(a[n].split(','))
        pts_list.append(b)

    for m in range(f.shape[0]):      # knot label
        a1 = f[m,0]
        knotid.append(a1)      
        
    for l in range(len(knotid)):    # full
        base_knots1[knotid[l]]=np.array([np.array(pts_list[l],dtype='d')*sf][0],dtype='d')
    
    del pts_list,f,knotid,b,a1,a
    gc.collect()
    
    return base_knots1