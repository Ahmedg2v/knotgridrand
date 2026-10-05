import numpy as np
import random
from .knot_structures import dictify, dictify_strand, array_knot, update_dict_add, update_dict_remove, update_dict_add_batch
from .grid import grid_path_3d
from .energy import energy_calc, energy_calc_local

# BFACF

def BFACF(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, theta):
    '''Run the BFACF algorithm.
    
    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    fwd : dict
        Maps each point index to the next point along the knot.

    bwd : dict
        Maps point indices to the previous point along the knot.

    gdc : dict
        Maps coordinate tuples to point indices.

    next_key : int
        New index not in the original dictionary.

    bx, by, bz, b_mat : array-like
        Cardinal vectors and 3D identity matrix, used for calculations.

    energy : float
        Energy of the initial configuration.

    T : float
        Temperature-like parameter for the probabilistic acceptance rules.

    alpha : float
        Length coefficient in energy calculation.

    beta : float
        Bends coefficient in energy calculation.

    gamma : float
        Distance to center of mass coefficient in energy calculation.

    theta : float
        Scaling factor for the probabilstic acceptance rules.
    
    Returns
    -------
    gd, fwd, bwd, gdc : dicts
        Updated knot dictionaries.

    next_key : int
        New key to be used as input for the next run.

    energy : float
        Energy of the new configuration.
    '''
    next_key = next_key+2
    idx = random.choice(list(gd.keys()))

    a1 = gd[idx]
    a2 = gd[fwd[idx]]
    a_1 = gd[bwd[idx]]
    a_2 = gd[bwd[bwd[idx]]]
    a3 = gd[fwd[fwd[idx]]]

    #b = ((a2-a1)+1)%2

    possible_moves = [[a1+bx,a2+bx],
                      [a1-bx,a2-bx],
                      [a1+by,a2+by],
                      [a1-by,a2-by],
                      [a1+bz,a2+bz],
                      [a1-bz,a2-bz]]

    die0 = random.choice(range(6))
    picked_move = possible_moves[die0]

    tup0 = f"{picked_move[0][0]},{picked_move[0][1]},{picked_move[0][2]}"
    tup1 = f"{picked_move[1][0]},{picked_move[1][1]},{picked_move[1][2]}"

    check0 = tup0 not in gdc
    check1 = tup1 not in gdc

    prev0 = f"{a_1[0]},{a_1[1]},{a_1[2]}" == tup0
    next1 = f"{a3[0]},{a3[1]},{a3[2]}" == tup1

    if check0 and check1: # no crossing
        delta_e = alpha*2 + beta*2 
        energy1 = energy+delta_e

        prob = np.exp(-theta*energy1/T)
        #flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0],1)
        
        flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0], 
                       np.where(delta_e < 0, random.choices([0,1],[prob,1-prob],k=1)[0], 1))
        
        if flip==1:
            gd, fwd, bwd, gdc = update_dict_add(next_key, idx, picked_move[0], gd, fwd, bwd, gdc)

            gd, fwd, bwd, gdc = update_dict_add(next_key+1, next_key, picked_move[1], gd, fwd, bwd, gdc)

            energy = energy1
        else:
            pass        

    elif prev0 and next1:
        delta_e = -alpha*2 - beta*2 
        energy1 = energy+delta_e

        prob = np.exp(-theta*energy1/T)
        #flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0],1)
        
        flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0], 
                       np.where(delta_e < 0, random.choices([0,1],[prob,1-prob],k=1)[0], 1))

        if flip==1:                
            gd, fwd, bwd, gdc = update_dict_remove(fwd[idx], gd, fwd, bwd, gdc)

            gd, fwd, bwd, gdc = update_dict_remove(idx, gd, fwd, bwd, gdc)

            energy = energy1
        else:
            pass

    elif prev0 and check1:
        local_pts_old = [a_2, a_1, a1, a2, a3]

        local_pts_new = [a_2, a_1, picked_move[1], a2, a3]  

        delta_e = energy_calc_local(local_pts_new, gd, alpha, beta, gamma)-energy_calc_local(local_pts_old, gd, alpha, beta, gamma)
        energy1 = energy+delta_e

        prob = np.exp(-theta*energy1/T)
        #flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0],1)
        
        flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0], 
                       np.where(delta_e < 0, random.choices([0,1],[prob,1-prob],k=1)[0], 1))

        if flip==1:
            gd, fwd, bwd, gdc = update_dict_add(next_key, bwd[idx], picked_move[1], gd, fwd, bwd, gdc)

            gd, fwd, bwd, gdc = update_dict_remove(idx, gd, fwd, bwd, gdc)

            energy = energy1
        else:
            pass

    elif next1 and check0:
        local_pts_old = [a_2, a_1, a1, a2, a3]

        local_pts_new = [a_2, a_1, picked_move[0], a2, a3]  

        delta_e = energy_calc_local(local_pts_new, gd, alpha, beta, gamma)-energy_calc_local(local_pts_old, gd, alpha, beta, gamma)
        energy1 = energy+delta_e

        prob = np.exp(-theta*energy1/T)
        #flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0],1)
        
        flip = np.where(delta_e > 0, random.choices([0,1],[1-prob,prob],k=1)[0], 
                       np.where(delta_e < 0, random.choices([0,1],[prob,1-prob],k=1)[0], 1))

        if flip==1:
            gd, fwd, bwd, gdc = update_dict_add(next_key, idx, picked_move[0], gd, fwd, bwd, gdc)

            gd, fwd, bwd, gdc = update_dict_remove(fwd[idx], gd, fwd, bwd, gdc)

            energy = energy1
        else:
            pass

    next_key+=1 
    
    return gd, fwd, bwd, gdc, next_key, energy

# LARGE SCALE SLIDE MOVES

def rng_slide(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, phi):
    '''Randomly pick a strand or random length along the knot and translate it a random 
    distance distance in a random cardinal direction with a probability to make a second translation. 
    
    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    fwd : dict
        Maps each point index to the next point along the knot.

    bwd : dict
        Maps point indices to the previous point along the knot.

    gdc : dict
        Maps coordinate tuples to point indices.

    next_key : int
        New index not in the original dictionary.

    bx, by, bz, b_mat : array-like
        Cardinal vectors and 3D identity matrix, used for calculations.

    energy : float
        Energy of the initial configuration.

    T : float
        Temperature-like parameter for the probabilistic acceptance rules.

    alpha : float
        Length coefficient in energy calculation.

    beta : float
        Bends coefficient in energy calculation.

    gamma : float
        Distance to center of mass coefficient in energy calculation.

    phi : float
        Scaling factor for the probabilstic acceptance rules.
    
    Returns
    -------
    gd, fwd, bwd, gdc : dicts
        Updated knot dictionaries.

    next_key : int
        New key to be used as input for the next run.

    energy : float
        Energy of the new configuration.
    '''
    worked=False
    next_key0 = next_key+5
    length=len(gd)
    die0 = random.sample(b_mat,2)
    axis = die0[0]*random.choice([1,-1])
    dist = random.choice(range(5,20)) 
    start_pt = random.choice(list(gd.keys()))    #idx of the 1st point afterwhich the strand starts
    #die1 = random.choice(range(int(0.01*length),int(0.1*length)))
    die1 = random.choice(range(5,30))

    strand_idx = []
    strand_idx.append(fwd[start_pt])

    for u in range(die1):
        strand_idx.append(fwd[strand_idx[u]])

    end_pt=fwd[strand_idx[-1]]
    strand_pts = []

    for l in range(die1):
        strand_pts.append(gd[strand_idx[l]])

    strand_pts=np.array(strand_pts)
    ghost_check = [f"{v[0]},{v[1]},{v[2]}" for v in strand_pts]
    check_str = []         

    for d1 in range(1,dist):
        current=strand_pts+(d1*axis)
        comp=[f"{v[0]},{v[1]},{v[2]}" for v in current]
        check_str.append(any(p in gdc and p not in ghost_check for p in comp))
        
    new_strand0=strand_pts+(dist*axis)
    
    if random.choice([True, False]):
        axis2 = die0[1]*random.choice([1,-1])
        dist2 = random.choice(range(5,20))

        for d2 in range(1,dist2):
            current=new_strand0+(d2*axis2)
            comp=[f"{v[0]},{v[1]},{v[2]}" for v in current]
            check_str.append(any(p in gdc and p not in ghost_check for p in comp))

        new_strand1=new_strand0+(dist2*axis2)

        leg1 = grid_path_3d(gd[strand_idx[0]], new_strand0[0])
        leg1_= grid_path_3d(leg1[-1], new_strand1[0])[1:]
        leg2_= grid_path_3d(new_strand1[-1],new_strand0[-1])
        leg2 = grid_path_3d(leg2_[-1], gd[strand_idx[-1]])[1:]
        new_strand = np.vstack((leg1, leg1_,new_strand1[1:-1], leg2_,leg2))

    else:
        leg1 = grid_path_3d(gd[strand_idx[0]], new_strand0[0])
        leg2 = grid_path_3d(new_strand0[-1], gd[strand_idx[-1]])
        new_strand = np.vstack((leg1, new_strand0[1:-1], leg2))

    check_str = np.array(check_str)
    new_keys = [f"{p[0]},{p[1]},{p[2]}" for p in new_strand]
    
    if len(set(new_keys)) != len(new_keys) or any(k in gdc for k in new_keys[1:-1]) or any(check_str==True):
        pass

    else:
        delta_e=energy_calc_local(new_strand, gd, alpha, beta, gamma)-energy_calc_local(strand_pts, gd, alpha, beta, gamma)
        energy1=energy+delta_e
        prob = np.exp(-phi*energy1/T)
        #print(prob)

        if random.choices([0,1],[1-prob,prob],k=1)[0]==1:
            sd, sf, sb, sdc = dictify_strand(new_strand[:-1], next_key0)
            #print(new_strand)
            for n0 in range(len(strand_idx)-1):
                gd, fwd, bwd, gdc = update_dict_remove(strand_idx[n0], gd, fwd, bwd, gdc)

            gd, fwd, bwd, gdc = update_dict_add_batch(start_pt, next_key0, next_key0+len(sd)-1, sd, sf, sb, sdc, gd, fwd, bwd, gdc)
            
            worked=True
            next_key=next_key0+len(new_strand)+10
            
        else:
            pass
            
    return gd, fwd, bwd, gdc, next_key, energy, worked

# THE RANDOMISATION LOOP

def randomise(gd, fwd, bwd, gdc, t=40000, next_key=10000, T=750, mod=4, alpha=6, beta=2, gamma=0, theta=1, phi=1, r_states=False, r_lengths=False, r_energies=False):
    '''Randomise a given knot using the BFACF algorithm and large scale slide moves.
    
    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    fwd : dict
        Maps each point index to the next point along the knot.

    bwd : dict
        Maps point indices to the previous point along the knot.

    gdc : dict
        Maps coordinate tuples to point indices.

    t : int
        Number of time steps.

    next_key : int
        New index not in the original dictionary.

    T : float
        Temperature-like parameter for the probabilistic acceptance rules.

    mod : int
        Number of time steps after which a slide move is attempted.

    alpha : float
        Length coefficient in energy calculation.

    beta : float
        Bends coefficient in energy calculation.

    gamma : float
        Distance to center of mass coefficient in energy calculation.

    theta, phi : float
        Scaling factors for the probabilstic acceptance rules for BFACF and Slide moves.
    
    Returns
    -------
    gd, fwd, bwd, gdc : dict
        Updated knot dictionaries.

    next_key : int
        New key to be used as input in case of running the randomisation again without starting over.

    energies, slides, lengths, states : list
        Optional lists that store the knot quantities/states as needed.
    '''
    energy = energy_calc(gd, fwd, bwd, alpha, beta, gamma)
    
    bx = np.array([1,0,0])
    by = np.array([0,1,0])
    bz = np.array([0,0,1])
    b_mat = [bx,by,bz]
    lengths=[]
    states=[]
    energies=[]

    for m in range(t):
        gd, fwd, bwd, gdc, next_key, energy = BFACF(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, theta)

        if m%mod==0:
            gd, fwd, bwd, gdc, next_key, energy, worked = rng_slide(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, phi)
            #if worked:
                #print(m)
        if r_states:
            states.append(array_knot(gd,fwd))
        if r_lengths:
            lengths.append(len(gd))
        if r_energies:
            energies.append(energy)

    return gd, fwd, bwd, gdc, next_key, states,energies, lengths