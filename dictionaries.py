import numpy as np

# DEFINING DICTIONARY LOGIC AND CONVERTING BETWEEN ARRAYS AND DICTS

def dictify(x):
    '''Convert ordered array of points into dictionaries
    
    Parameters:
    x: array pts
    
    Returns:
    grid_dict: dictionary containing index:coordinate pairs     
    forward_map: dictionary storing order of pts 
    backward_map: dictionary storing order of pts in the opposite direction to forward_map
    grid_dict_check: dictionary containing coordinate:index pairs
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
    '''Convert ordered array of points into dictionaries
    
    Parameters:
    x: array pts
    key: starting index of strand pts (must not be in the original knot)
    
    Returns:
    grid_dict: dictionary containing index:coordinate pairs     
    forward_map: dictionary storing order of pts 
    backward_map: dictionary storing order of pts in the opposite direction to forward_map
    grid_dict_check: dictionary containing coordinate:index pairs
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
    '''Convert dictionary representation into ordered array
    
    Parameters:
    gd: dictionary containing index:coordinate pairs
    fwd: forward connectivity map
    
    Returns:
    grid_points: numpy array
    '''
    tlist=[]
    goto=list(fwd.keys())[0]
    for i in gd.keys():
        tlist.append(gd[fwd[goto]])
        goto=fwd[goto]

    grid_points=np.array(tlist)
    
    return grid_points

# UPDATING DICTIONARIES

def update_dict_add(x, index, move, grid_dict, cache_dict, cache_back, grid_dict_check):
    '''add a new point after the point at index'''
    
    cache_back[cache_dict[index]]=x
    cache_dict[x]=cache_dict[index]
    cache_dict[index]=x
    cache_back[x]=index
    grid_dict[x]=move
    grid_dict_check[f"{move[0]},{move[1]},{move[2]}"]=x
    
    return grid_dict, cache_dict, cache_back, grid_dict_check

def update_dict_remove(index, grid_dict, cache_dict, cache_back, grid_dict_check):
    '''remove the point at index'''
    
    cache_dict[cache_back[index]]=cache_dict[index]
    cache_back[cache_dict[index]]=cache_back[index]
    cache_back.pop(index)
    cache_dict.pop(index)
    gpop = grid_dict[index]
    grid_dict_check.pop(f"{gpop[0]},{gpop[1]},{gpop[2]}")
    grid_dict.pop(index)
    
    return grid_dict, cache_dict, cache_back, grid_dict_check

def update_dict_add_batch(index, str_start, str_end, str_dict, str_cache, str_back, str_check, grid_dict, cache_dict, cache_back, grid_dict_check):
    '''stitch a strand dict to the knot dict at index'''
    
    cache_dict.update(str_cache)
    cache_back.update(str_back)
    cache_back[cache_dict[index]]=str_end
    cache_dict[str_end]=cache_dict[index]
    cache_dict[index]=str_start
    cache_back[str_start]=index
    grid_dict.update(str_dict)
    grid_dict_check.update(str_check)
    
    return grid_dict, cache_dict, cache_back, grid_dict_check