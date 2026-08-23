import numpy as np

# ENERGY CALCULATION

def energy_calc(gd, fwd, bwd, alpha, beta, gamma):
    '''calculate the energy of a certain configuration.

    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    fwd : dict
        Maps each point index to the next point along the knot.

    bwd : dict
        Maps point indices to the previous point along the knot.

    alpha : float
        Length coefficient.

    beta : float
        Bends coefficient.

    gamma : float
        Distance to center of mass coefficient.

    Returns
    -------
    energy : float
        Energy of the configuration.
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
    '''calculate the energy of a section of the knot.

    Parameters
    ----------
    x : int
        Index of the point around which the local energy is calculated.

    gd : dict
        Knot dictionary.

    alpha : float
        Length coefficient.

    beta : float
        Bends coefficient.

    gamma : float
        Distance to center of mass coefficient.

    window : int, default=5
        Number of points preceding and proceeding x used in the energy calculation.

    Returns
    -------
    energy : float
        Energy of the local section.
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