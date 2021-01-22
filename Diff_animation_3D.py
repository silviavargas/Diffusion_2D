# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 11:48:48 2020
@author: Silvia Vargas
"""
import numpy as np
print('numpy: '+np.version.full_version)
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
print('matplotlib: '+matplotlib.__version__)
from Diffusion_2D import Diffusion
import sys
from sys import argv
import configparser

#%% 
              
config = configparser.ConfigParser()
config.read(sys.argv[1])

dx = config.get('parameters', 'dx')
nu = config.get('parameters', 'nu')
kind = config.get('parameters', 'kind')
nt = config.get ('parameters', 'nt')
L = config.get ('parameters', 'L')

name = config.get('gif names', 'name')

dx = float(dx)
nu = float(nu)
nt = int(nt)
L = float(L)

Diff1=Diffusion(dx, nu, kind, nt, L)
#%%
fps = 10 # frame per sec
Diff1.evolve_ts()

print ('Total time:' , round(Diff1.dt*nt*1000), 'ms')

x = np.linspace(0,Diff1.L, Diff1.nx)
y = np.linspace(0,Diff1.L, Diff1.ny)
X, Y = np.meshgrid(x,y)

    
def update_plot(frame_number, zarray, plot):
     plot[0].remove()
     plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_number], cmap="plasma")
     ax.set_title("Time {} [ms]".format(round(Diff1.dt*frame_number*1000)), fontsize=12)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plot = [ax.plot_surface(X, Y, Diff1.u[:,:,0], color='0.75', rstride=1, cstride=1)]
ax.set_zlim(0,1.1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ani = animation.FuncAnimation(fig, update_plot, Diff1.nt, fargs=(Diff1.u, plot), interval=1000/fps)

fn = name
ani.save(fn+'.gif',writer='imagemagick',fps=fps)