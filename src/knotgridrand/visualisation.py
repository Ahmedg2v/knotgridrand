import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import draw,pause,close
from mpl_toolkits.mplot3d import Axes3D
from .knot_structures import array_knot
import matplotlib.animation as animation

def plot_knot(gd, cd):
    '''Plot a knot in 3D.
    
    Parameters
    ----------
    gd : dict
        Maps point indices to coordinates.

    cd : dict
        Maps each point index to the next point along the knot.
    '''
    grid_points = array_knot(gd, cd)
    #line plot

    fig = plt.figure(figsize=[8,6])
    ax = plt.axes(projection='3d')
    #ax.grid(False)
    
    #ax.set_axis_off()
    #ax.grid(True)
    plot_points=np.vstack((grid_points,grid_points[0]))
    ax.plot3D(*plot_points.T, 'green')
    #plt.axis('off')
    plt.show()
    
def plot_knot_array(gp):
    '''Plot a knot in 3D.
    
    Parameters
    ----------
    gp : array-like, shape (N,3)
        Ordered coordinates describing the knot.
    '''
    grid_points = gp
    #line plot

    fig = plt.figure(figsize=[8,6])
    ax = plt.axes(projection='3d')
    
    plot_points=np.vstack((grid_points,grid_points[0]))
    ax.plot3D(*plot_points.T, 'green')
    #plt.xlim(-2, 20)
    #plt.ylim(-20, 15)
    #plt.axis('off')
    plt.show()
    
def graph_avg(x,label):
    '''Plot a graph of the average of a list of values.
    
    Parameters
    ----------
    x : list
        List of values.

    label : str
        Label for the y-axis.
    '''
    data=np.array(x)
    t = np.arange(0,len(x))
    fig = plt.figure(figsize=(12,6))
    
    plt.plot(t, data, 'g')

    plt.xlabel('time step', fontsize=14)
    plt.ylabel(label, fontsize=14)

    plt.show()

def animate_knot_evolution(states, coord_min, coord_max, step=10, interval=1):
    '''
    Animates a list of (N, 3) knot states using pre-known axes limits.

    Parameters
    ----------
    states : list of arrays, shape (N, 3)
        List of knot states over time.

    coord_min : list of float
        Minimum coordinates for each axis.

    coord_max : list of float
        Maximum coordinates for each axis.

    step : int
        Frame skip (stride) for animation.

    interval : int
        Delay between frames in milliseconds.
    '''
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Apply your known limits directly
    ax.set_xlim([coord_min[0], coord_max[0]])
    ax.set_ylim([coord_min[1], coord_max[1]])
    ax.set_zlim([coord_min[2], coord_max[2]])
    
    ax.set_title("Knot Evolution")
    ax.axis('on') 
    
    line, = ax.plot([], [], [], lw=1, color='green')
    
    # Apply the frame skip (stride)
    plot_states = states[1::step]
    #print(f"Animating {len(plot_states)} frames (Step size: {step})...")

    # The update function
    def update(frame_index):
        current_state = plot_states[frame_index]
        
        # Slicing the (N, 3) array
        x = current_state[:, 0]
        y = current_state[:, 1]
        z = current_state[:, 2]
        
        line.set_data(x, y)
        line.set_3d_properties(z)
        
        return [line]

    ani = animation.FuncAnimation(
        fig, 
        update, 
        frames=len(plot_states), 
        interval=interval, 
        blit=False, 
        repeat=True
    )
    
    plt.show()
