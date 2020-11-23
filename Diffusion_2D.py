# -*- coding: utf-8 -*-
"""
author: Silvia Vargas

"""
import numpy as np

def diffuse(I, nu, L, F, T): #tolto dt come parametro
    import time;  t0 = time.clock()  # For measuring the CPU time
    dx,dy= L* 0.01
    dx2, dy2 = dx*dx, dy*dy
    dt = F * (dx2 * dy2) / (nu * (dx2 + dy2))
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    Nx = int(round(L/dx))
    Ny = int(round(L/dy))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space 
    y = np.linspace(0, L, Ny+1)        
    
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = t[1] - t[0]
    
    u   = np.zeros((Nx+1, Ny+1, Nt+1))      # solution array     
     
# Set initial condition
    for i in range(0,Nx+1):
     for j in range (0,Ny+1):
         u[i,j,0] = I[x[i],y[j]]#modificato da funz a array
 
#Run through Nt timesteps    
    for n in range(1, Nt):       
      u[1:-1, 1:-1, n] = u[1:-1, 1:-1, n-1] + nu * dt * (
          (u[2:, 1:-1, n-1] - 2*u[1:-1, 1:-1, n-1] + u[:-2, 1:-1,n-1])/dx2
          + (u[1:-1, 2: ,n-1] - 2*u[1:-1, 1:-1, n-1] + u[1:-1, :-2, n-1])/dy2 )
      
    
    t1 = time.clock()
    return u, t1-t0 #u_n? 



