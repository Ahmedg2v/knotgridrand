import numpy as np
import random
import pyknotid.spacecurves as pks
from pyknotid.spacecurves import Knot
import pyknotid.invariants as pki
import pyknotid as pk
from scipy.spatial.transform import Rotation as R
from dictionaries import array_knot

def sta_writhe(grid_dict,forward_map,width):
    with np.errstate(divide='ignore', invalid='ignore'):
        wknot = array_knot(grid_dict,forward_map)
        wtang = np.vstack((wknot[1:],wknot[0]))-wknot
        dots=[]
        mods=[]
        for kl in range(len(wknot)):
            vc = wknot[kl]-wknot
            cr = np.cross(wtang[kl], wtang)
            print(vc.shape,cr.shape)
            
            dot=[]
            for m in range(len(wknot)):
                dot.append(np.dot(cr[m], vc[m]))
                print(dot[-1],dot)
            dots.append(dot)
            mod = np.linalg.norm(vc,axis=1)**3
            mods.append(mod)

        dots=np.array(dots)
        mods=np.array(mods)
        frac=dots/mods
        frac[np.isnan(frac)]=0
        wr = np.sum(frac,axis=1).reshape(-1, 1)

        offset=30+width
        wr_padded = np.vstack((wr[-offset:],wr,wr[:offset]))
        sta=[]
        for i in range(len(wknot)):
            lim1, lim2 = int(offset+i-width/2), int(offset+i+width/2)
            sta.append(np.sum(wr_padded[lim1:lim2]))
        sta=np.array(sta).reshape(-1, 1)
    
    return sta

def rotate(pts):
    quat = R.random().as_quat()    # (x, y, z, w)
    rot = R.from_quat(quat)
    centre = np.array([np.average(pts[:,0]),np.average(pts[:,1]),np.average(pts[:,2])])
    return rot.apply(pts-centre)+centre
    
def alex(pts,var):
    quat = R.random().as_quat()    # (x, y, z, w)
    rot = R.from_quat(quat)
    rotated_pts = rot.apply(pts)
    K = Knot(rotated_pts)
    ktry=K.gauss_code()
    alex=pki.alexander(ktry, variable=var, simplify=True, mode='python')
    
    return alex

def add_pts(arr, target):
    #arr=array_knot(grid_dict,forward_map)          #not for bluebear - rand knots already converted to arrays 
    mids=(np.vstack((arr[1:],arr[0]))+arr)/2

    mask=np.full((len(mids),1),np.nan)
    mask[0:target-len(arr)]=1
    p=np.random.permutation(len(mask))
    mask=mask[p]
    mids=mids*mask

    tmp=np.zeros((len(mids)+len(arr),3))
    tmp[::2]=tmp[::2]+arr
    tmp[1::2]=tmp[1::2]+mids

    new_arr = tmp[np.isfinite(tmp).all(axis=1)]

    if (len(new_arr)<target):
        new_arr = add_pts(new_arr, target)
    else:
        new_arr = new_arr 

    return new_arr

def haus_like_area(hknot1,hknot2):
    hdist1=[]
    hdist2=[]

    for h1 in range(len(hknot1)):
        hdist1.append(np.min(np.linalg.norm(hknot1[h1]-hknot2,axis=1)))

    for h2 in range(len(hknot2)):
        hdist2.append(np.min(np.linalg.norm(hknot2[h2]-hknot1,axis=1)))

    hdist1=np.array(hdist1)
    hdist2=np.array(hdist2)
    area1 = np.sum(hdist1)
    area2 = np.sum(hdist2)
    
    return np.max([area1,area2])

def haus(hknot1,hknot2):
    hdist1=[]
    hdist2=[]

    for h1 in range(len(hknot1)):
        hdist1.append(np.min(np.linalg.norm(hknot1[h1]-hknot2,axis=1)))

    for h2 in range(len(hknot2)):
        hdist2.append(np.min(np.linalg.norm(hknot2[h2]-hknot1,axis=1)))

    hdist1=np.array(hdist1)
    hdist2=np.array(hdist2)
    
    return np.max([np.max(hdist1),np.max(hdist2)])