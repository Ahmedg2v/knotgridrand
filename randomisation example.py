import numpy as np
from dictionaries import dictify
from rasterise import get_grid
from randomise import randomise
from graphs import plot_knot
#base_knots = generate_base_knots_dict('3d-coordinates.csv')
#np.savez_compressed("base_knots_10.npz", **base_knots)
# t, next_key, T, alpha, beta, gamma, theta, phi
base_knots = np.load("base_knots_20.npz")
grid_points = get_grid('3_1', base_knots)
grid_dict, forward_map, backward_map, grid_dict_check = dictify(grid_points)
grid_dict, forward_map, backward_map, grid_dict_check, next_key, energies, slides, lengths, states= randomise(grid_dict, forward_map, backward_map, grid_dict_check)
plot_knot(grid_dict,forward_map)