import numpy as np
import time
from src.knotgridrand.datasets import generate_base_knots_dict
from src.knotgridrand.knot_structures import dictify
from src.knotgridrand.grid import get_grid
from src.knotgridrand.randomisation import randomise
from src.knotgridrand.visualisation import plot_knot
from src.knotgridrand.visualisation import plot_knot_array
from src.knotgridrand.knot_structures import array_knot
from src.knotgridrand.visualisation import animate_knot_evolution
#base_knots = generate_base_knots_dict('3d-coordinates.csv')
#np.savez_compressed("base_knots_10.npz", **base_knots)
base_knots = np.load("base_knots_sf20_10crossings.npz")
plot_knot_array(base_knots["3_1"])
grid_points = get_grid("3_1", base_knots)

plot_knot_array(grid_points)
grid_dict, forward_map, backward_map, grid_dict_check = dictify(grid_points)

grid_dict, forward_map, backward_map, grid_dict_check, states, energies, lengths= randomise(grid_dict, forward_map, backward_map, grid_dict_check, T=1000, r_step=20, strand_length=30, r_states=True)
plot_knot(grid_dict,forward_map)

# ANIMATION EXAMPLE USING THE STATES RETURNED FROM RANDOMISATION OF A TREFOIL KNOT
# keep in mind the effective no. of states between consecutive frames is r_step*step
animate_knot_evolution(states, coord_min=None, coord_max=None, step=1, interval=1)

#np.savez_compressed("states1.npz", *states, allow_pickle=True)
#data = np.load("states.npz", allow_pickle=True)
#states = [data[key] for key in data.files]
