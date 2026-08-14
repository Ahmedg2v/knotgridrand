import numpy as np

# DEFINING DICTIONARY LOGIC AND CONVERTING BETWEEN ARRAYS AND DICTS

def dictify(x):
    '''Convert ordered array of points into dictionaries.
    
    Parameters
    ----------
    x : array-like, shape (N,3)
        Ordered coordinates describing the knot.

    Returns
    -------
    grid_dict: dict
        Maps point indices to coordinates.

    forward_map : dict
        Maps each point index to the next point along the knot.

    backward_map : dict
        Maps each point index to the previous point along the knot.

    grid_dict_check : dict
        Maps coordinate tuples to point indices.
    '''
    grid_dict = {}
    forward_map = {}
    backward_map = {}

    for j in range(len(x)-1):
        grid_dict[j] = np.array([x[j][0],x[j][1],x[j][2]])
        forward_map[j] = (j+1)%(len(x)-1)
        backward_map[j] = (j-1)%(len(x)-1)
        
    grid_dict_check = { f"{v[0]},{v[1]},{v[2]}": k for k,v in grid_dict.items() }

    return grid_dict, forward_map, backward_map, grid_dict_check

def dictify_strand(x, key):
    '''Convert ordered array of points into dictionaries.
    
    Parameters
    ----------
    x : array-like, shape (N,3)
        Ordered coordinates describing the strand.

    key : int
        Starting index of strand points (must not be in the original knot).

    Returns
    -------
    grid_dict : dict
        Maps point indices to coordinates.

    forward_map : dict
        Maps each point index to the next point along the strand.

    backward_map : dict
        Maps each point index to the previous point along the strand.

    grid_dict_check : dict
        Maps coordinate tuples to point indices.
    '''
    grid_dict = {}
    forward_map = {}
    backward_map = {}

    for j in range(len(x)):
        grid_dict[j+key] = np.array([x[j][0],x[j][1],x[j][2]])
        forward_map[j+key] = (j+1)%len(x)+key
        backward_map[j+key] = (j-1)%len(x)+key
        
    grid_dict_check = { f"{v[0]},{v[1]},{v[2]}": k for k,v in grid_dict.items() }
    return grid_dict, forward_map, backward_map, grid_dict_check

def array_knot(gd, fwd):
    '''Convert dictionary representation into ordered array.
    
    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    fwd : dict
        Maps each point index to the next point along the knot.

    Returns
    -------
    grid_points : array-like, shape (N,3)
        Ordered coordinates describing the knot.
    '''
    tlist=[]
    goto=list(fwd.keys())[0]
    for i in gd.keys():
        tlist.append(gd[fwd[goto]])
        goto=fwd[goto]

    grid_points=np.array(tlist)
    
    return grid_points

# UPDATING DICTIONARIES

def update_dict_add(x, index, move, grid_dict, forward_map, backward_map, grid_dict_check):
    '''add a new point after the point at index.
    
    Parameters
    ----------
    x : int
        Index of the new point to be added.

    index : int
        Index of the point after which to add the new point.

    move : array-like, shape (3,)
        Coordinates of the new point.

    grid_dict : dict
        Maps point indices to coordinates.

    forward_map : dict
        Maps each point index to the next point along the knot.

    backward_map : dict
        Maps each point index to the previous point along the knot.

    grid_dict_check : dict
        Maps coordinate tuples to point indices.

    Returns
    -------
    grid_dict : dict
        Updated point indices to coordinates mapping.

    forward_map : dict
        Updated forward connectivity map.

    backward_map : dict
        Updated backward connectivity map.

    grid_dict_check : dict
        Updated coordinate tuples to point indices mapping.
    '''
    
    backward_map[forward_map[index]]=x
    forward_map[x]=forward_map[index]
    forward_map[index]=x
    backward_map[x]=index
    grid_dict[x]=move
    grid_dict_check[f"{move[0]},{move[1]},{move[2]}"]=x
    
    return grid_dict, forward_map, backward_map, grid_dict_check

def update_dict_remove(index, grid_dict, forward_map, backward_map, grid_dict_check):
    '''remove the point at index.
    
    Parameters
    ----------
    index : int
        Index of the point to be removed.

    grid_dict : dict
        Maps point indices to coordinates.

    forward_map : dict
        Maps each point index to the next point along the knot.

    backward_map : dict
        Maps each point index to the previous point along the knot.

    grid_dict_check : dict
        Maps coordinate tuples to point indices.

    Returns
    -------
    grid_dict : dict
        Updated point indices to coordinates mapping.

    forward_map : dict
        Updated forward connectivity map.

    backward_map : dict
        Updated backward connectivity map.

    grid_dict_check : dict
        Updated coordinate tuples to point indices mapping.
    '''
    
    forward_map[backward_map[index]]=forward_map[index]
    backward_map[forward_map[index]]=backward_map[index]
    backward_map.pop(index)
    forward_map.pop(index)
    gpop = grid_dict[index]
    grid_dict_check.pop(f"{gpop[0]},{gpop[1]},{gpop[2]}")
    grid_dict.pop(index)
    
    return grid_dict, forward_map, backward_map, grid_dict_check

def update_dict_add_batch(index, str_start, str_end, str_dict, str_cache, str_back, str_check, grid_dict, forward_map, backward_map, grid_dict_check):
    '''stitch a strand to the knot at index.

    Parameters
    ----------
    index : int
        Index of the point to which the strand will be attached.
    
    str_start : int
        Index of the first point in the strand.

    str_end : int
        Index of the last point in the strand.

    str_dict : dict
        Maps point indices to coordinates for the strand.

    str_cache : dict
        Maps each point index to the next point along the strand.

    str_back : dict
        Maps each point index to the previous point along the strand.

    str_check : dict
        Maps coordinate tuples to point indices for the strand.

    grid_dict: dict
        Maps point indices to coordinates for the knot.
    
    forward_map : dict
        Maps each point index to the next point along the knot.

    backward_map : dict
        Maps each point index to the previous point along the knot.

    grid_dict_check : dict
        Maps coordinate tuples to point indices for the knot.

    Returns
    -------
    grid_dict: dict
        Maps point indices to coordinates.

    forward_map : dict
        Maps each point index to the next point along the knot.

    backward_map : dict
        Maps each point index to the previous point along the knot.

    grid_dict_check : dict
        Maps coordinate tuples to point indices.
    '''
    
    forward_map.update(str_cache)
    backward_map.update(str_back)
    backward_map[forward_map[index]]=str_end
    forward_map[str_end]=forward_map[index]
    forward_map[index]=str_start
    backward_map[str_start]=index
    grid_dict.update(str_dict)
    grid_dict_check.update(str_check)
    
    return grid_dict, forward_map, backward_map, grid_dict_check