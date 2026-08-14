import numpy as np
import random

# PROJECTING BASE KNOTS ONTO THE 3D LATTICE

def get_grid(x, base_knots):
    '''Generate a rasterised version of a specific knot.
    
    Parameters
    ----------
    x: string
        ID of the knot
    base_knots: dict
        Dictionary containing the database knots

    Returns
    -------
    grid_points: numpy array of shape (N,3) with integer values
        Rasterised version of the knot
    '''
    knot = base_knots[x]
    grid_points = []
    for i in range(len(knot)):
        path = grid_path_3d(knot[i], knot[(i+1)%len(knot)])
        grid_points.extend(path)

    grid_points=np.array(grid_points,dtype=int)    
    mask = np.any(np.diff(grid_points, axis=0) != 0, axis=1)
    # True to keep the first point
    mask = np.r_[True, mask]

    # result
    grid_points = grid_points[mask]
    
    grid_points=np.array(grid_points,dtype=float)
    
    b=np.unique(grid_points[:-1], axis=0, return_counts=True)     #count instances of same coordinate (b[1][k1,2]=2 at both pts)
    c=np.argwhere(b[1]-1==1)                                      #idx of duplicates in b

    for k in c:
        idx = np.argwhere(np.all(grid_points==b[0][k],axis=1))
        lim1, lim2 = int(idx[0]+1), int(idx[1]+1)
        grid_points[lim1:lim2]=np.nan
    
    grid_points = grid_points[np.isfinite(grid_points).all(axis=1)]
            
    grid_points = np.array(grid_points).astype(np.int64)
    
    return grid_points


def grid_path_3d(p0, p1):
    '''Rasterize the line segment connecting two points such that there are no diagonal steps.

    Parameters
    ----------
    p0: array-like, shape (3,)
        Start point.

    p1: array-like, shape (3,)
        End point.

    Returns
    -------
    Rasterised line segment as a numpy array of shape (N,3) with integer values
    '''
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)

    x, y, z = np.floor(p0).astype(int)
    x1, y1, z1 = np.floor(p1).astype(int)

    dx, dy, dz = p1 - p0
    dr = [dx,dy,dz]
    step = np.sign(dr).astype(int)

    inv_dr = []
    paradist = []
    for i in range(3):
        if dr[i] != 0:
            inv_dr.append(abs(1/dr[i]))
        else:
            inv_dr.append(np.inf)
    
    for j in range(3):
        if dr[j] == 0:
            paradist.append(np.inf)
        else:
            if step[j] > 0:
                dist = (np.floor(p0[j]) + 1 - p0[j])
            else:
                dist = (p0[j] - np.floor(p0[j]))

            paradist.append(dist * inv_dr[j])
        
    path = [(x, y, z)]
    
    nx = abs(x1 - x)
    ny = abs(y1 - y)
    nz = abs(z1 - z)

    while nx > 0 or ny > 0 or nz > 0:
        if paradist[0] <= paradist[1] and paradist[0] <= paradist[2]:
            paradist[0] += inv_dr[0]
            x += step[0]
            nx -= 1
        elif paradist[1] <= paradist[2]:
            paradist[1] += inv_dr[1]
            y += step[1]
            ny -= 1
        else:
            paradist[2] += inv_dr[2]
            z += step[2]
            nz -= 1

        path.append((x, y, z))

    return np.array(path, dtype=int)

