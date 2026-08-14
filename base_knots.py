import numpy as np
import pandas as pd
import gc

# IMPORTING A SUBSET OF BASE KNOTS

def import_base_knots(x):
    '''import a csv file containing the base knots and convert it into a numpy array.'''
    p05 = pd.read_csv(x)
    data_g = np.array(p05)

    return data_g

def generate_base_knots_dict(x, n=0, N=250, sf=20):
    '''generate a dictionary containing the specified range of knots.
    
    Parameters
    ----------
    x: csv file storing 3D coordinates of representative knots.

    n, N: int
        Range of the knots to be imported.

    sf: int
        Scaling factor multiplied by the coordinates to avoid float point inaccuracies during rasterisation.

    Returns
    -------
    base_knots: dict
        Dictionary containing 'knot id: coords array' pairs.
    '''
    # max N is 12965 in the knotinfo dataset
    f = import_base_knots(x)[n:N]
    pts_list=[]    
    knotid=[]
    base_knots = {}
    
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
        base_knots[knotid[l]]=np.array([np.array(pts_list[l],dtype='d')*sf][0],dtype='d')
    
    del pts_list,f,knotid,b,a1,a
    gc.collect()
    
    return base_knots