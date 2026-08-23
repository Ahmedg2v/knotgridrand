import numpy as np
from src.knotgridrand.datasets import generate_base_knots_dict
from src.knotgridrand.knot_structures import dictify
from src.knotgridrand.grid import get_grid
from src.knotgridrand.randomisation import randomise
from src.knotgridrand.visualisation import plot_knot
#base_knots = generate_base_knots_dict('3d-coordinates.csv')
#np.savez_compressed("base_knots_10.npz", **base_knots)
base_knots = np.load("base_knots_sf20_10crossings.npz")
grid_points = get_grid("3_1", base_knots)
grid_dict, forward_map, backward_map, grid_dict_check = dictify(grid_points)
grid_dict, forward_map, backward_map, grid_dict_check, next_key, energies, slides, lengths, states= randomise(grid_dict, forward_map, backward_map, grid_dict_check)
plot_knot(grid_dict,forward_map)