import numpy as np
from src.knotgridrand.datasets import generate_base_knots_dict
from src.knotgridrand.knot_structures import dictify
from src.knotgridrand.grid import get_grid
from src.knotgridrand.randomisation import randomise
from src.knotgridrand.visualisation import plot_knot
from src.knotgridrand.visualisation import plot_knot_array
from src.knotgridrand.knot_structures import array_knot
from src.knotgridrand.visualisation import animate_knot_evolution
base_knots = generate_base_knots_dict('3d-coordinates.csv')
#np.savez_compressed("base_knots_10.npz", **base_knots)
#base_knots = np.load("base_knots_sf20_10crossings.npz")
grid_points = get_grid("3_1", base_knots)
grid_dict, forward_map, backward_map, grid_dict_check = dictify(grid_points)
plot_knot_array(base_knots["3_1"])
grid_dict, forward_map, backward_map, grid_dict_check, next_key, states, energies, lengths= randomise(grid_dict, forward_map, backward_map, grid_dict_check, T=1000)
plot_knot(grid_dict,forward_map)

# ANIMATION EXAMPLE USING THE STATES RETURNED FROM RANDOMISATION OF A TREFOIL KNOT

#np.savez_compressed("states.npz", *states, allow_pickle=True)
data = np.load("states.npz", allow_pickle=True)
states = [data[key] for key in data.files]
animate_knot_evolution(states, coord_min=[-10, -50, -30], coord_max=[30, 0, 20], step=20, interval=1)
