import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.pyplot import draw,pause,close
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from dictionaries import array_knot

def plot_try(x):
    fig = plt.figure(figsize=[8,6])
    ax = plt.axes(projection='3d')
    mm = base_knots[x]
    x=[]
    y=[]
    z=[]

    for k in np.arange(0, len(mm)+1,1):
        x.append(mm[k%len(mm)][0])
        y.append(mm[k%len(mm)][1])
        z.append(mm[k%len(mm)][2])

    ax.plot3D(x, y, z, 'green')

    plt.show()

def plot_knot(gd, cd):
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
    data=np.array(x)
    t = np.arange(0,len(x))
    fig = plt.figure(figsize=(12,6))
    
    plt.plot(t, data, 'g')

    plt.xlabel('time step', fontsize=14)
    plt.ylabel(label, fontsize=14)

    plt.show()