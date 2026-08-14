import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import draw,pause,close
from mpl_toolkits.mplot3d import Axes3D
from dictionaries import array_knot

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