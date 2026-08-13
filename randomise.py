import numpy as np
import random
from dictionaries import dictify, dictify_strand, array_knot, update_dict_add, update_dict_remove, update_dict_add_batch
from rasterise import grid_path_3d

# ENERGY CALCULATION

def energy_calc(gd, fwd, bwd, alpha, beta, gamma):
    '''calculate the energy of a certain configuration
    
    Parameters:
    gd, fwd, bwd: knot dictionaries
    
    coefficients of different energy terms:
    alpha: length coeff
    beta: bends coeff
    gamma: distance to com coeff
    
    Returns: energy as a float
    '''
    dot_sum=0
    potential=0
    pts = np.array(list(gd.values()))
    
    x_mean = np.average(pts[:,0])
    y_mean = np.average(pts[:,1])
    z_mean = np.average(pts[:,2])
    
    for p in list(gd.keys()):
        bwdd = bwd[p]
        fwdd = fwd[p]
        
        gd0 = gd[bwdd]
        gd1 = gd[p]
        gd2 = gd[fwdd]
        
        potential+=np.sqrt((gd1[0]-x_mean)**2+
                           (gd1[1]-y_mean)**2+
                           (gd1[2]-z_mean)**2)
    
        dv0 = gd2-gd1
        dv1 = gd1-gd0

        dot_sum += np.dot(dv0, dv1)
    
    energy = alpha*len(gd)+beta*(len(gd)-dot_sum)+gamma*potential
    return energy

def energy_calc_local(x, gd, alpha, beta, gamma, window=5):
    '''calculate the energy of a section of the knot
    
    Parameters:
    x: index of the point around which the local energy is calculated
    window: number of pts preceding and proceeding x used in the energy calculation
    gd: knot dictionary
    
    coefficients of different energy terms:
    alpha: length coeff
    beta: bends coeff
    gamma: distance to com coeff
    
    Returns: energy as a float
    '''
    dot_sum=0
    potential=0
    pts = np.array(list(gd.values()))
    
    x_mean = np.average(pts[:,0])
    y_mean = np.average(pts[:,1])
    z_mean = np.average(pts[:,2])
    
    for p in range(window):  
        x0 = x[(p-1)%len(x)]
        x1 = x[p]
        x2 = x[(p+1)%len(x)]
        
        potential+=np.sqrt((x1[0]-x_mean)**2+
                           (x1[1]-y_mean)**2+
                           (x1[2]-z_mean)**2)
    
        dv0 = x2-x1
        dv1 = x1-x0

        dot_sum += np.dot(dv0, dv1)
    
    energy = alpha*len(x)+beta*(len(x)-dot_sum)+gamma*potential
    return energy

# BFACF

def BFACF(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, theta):
    '''Run the BFACF algorithm 
    
    Parameters:
    gd, fwd, bwd, gdc: knot dictionaries
    next_key: new index not in the original dict
    bx, by, bz, b_mat: cardinal vectors and 3D identity matrix, used for calculations
    energy: energy of the initial configuration
    T: temperature like parameter
    alpha, beta, gamma: energy coefficients
    theta: scaling factor for the probabilstic acceptance rules
    
    Returns:
    gd, fwd, bwd, gdc: updated knot dictionaries
    next_key: new key to be used as input for the next run
    energy: energy of the new configuration
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
    distance distance in a random cardinal direction with a probability to make a second translation 
    
    Parameters:
    gd, fwd, bwd, gdc: knot dictionaries
    next_key: new index not in the original dict
    bx, by, bz, b_mat: cardinal vectors and 3D identity matrix, used for calculations
    energy: energy of the initial configuration
    T: temperature like parameter
    alpha, beta, gamma: energy coefficients
    phi: scaling factor for the probabilstic acceptance rules
    
    Returns:
    gd, fwd, bwd, gdc: updated knot dictionaries
    next_key: new key to be used as input for the next run
    energy: energy of the new configuration
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

def randomise(gd, fwd, bwd, gdc, t=40000, next_key=10000, T=750, mod=4, alpha=6, beta=2, gamma=0, theta=1, phi=1):
    '''Randomise a given knot
    
    Parameters:
    gd, fwd, bwd, gdc: dictionaries storing the knot
    t: number of time steps
    next_key: new index not in the original dict
    T: temperature like parameter
    mod: number of time steps after which a slide move is attempted
    alpha, beta, gamma: energy parameters
    theta, phi: scaling factors for the probabilstic acceptance rules
    
    Returns:
    gd, fwd, bwd, gdc: updated knot dictionaries
    next_key: new key to be used as input in case of running the randomisation again without starting over
    energies, slides, lengths, states: optional lists that store the knot quantities/states as needed
    '''
    energy = energy_calc(gd, fwd, bwd, alpha, beta, gamma)
    
    bx = np.array([1,0,0])
    by = np.array([0,1,0])
    bz = np.array([0,0,1])
    b_mat = [bx,by,bz]
    lengths=[]
    states=[]
    energies=[]
    slides=[]
    for m in range(t):
        gd, fwd, bwd, gdc, next_key, energy = BFACF(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, theta)

        if m%mod==0:
            gd, fwd, bwd, gdc, next_key, energy, worked = rng_slide(gd, fwd, bwd, gdc, next_key, bx, by, bz, b_mat, energy, T, alpha, beta, gamma, phi)
            #if worked:
                #print(m)
        #states.append(array_knot(gd,fwd))
        #lengths.append(len(gd))
        #energies.append(energy)
    return gd, fwd, bwd, gdc, next_key, energies, slides, lengths, states