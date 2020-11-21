# -*- coding: utf-8 -*-
"""
author: Silvia Vargas

"""
import numpy as np

def diffuse(I, nu, dx, dy, L, dt, F, T): 
    import time;  t0 = time.clock()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    Nx = int(round(L/dx))
    Ny = int(round(L/dy))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space 
    y = np.linspace(0, L, Ny+1)        
    dx2, dy2 = dx*dx, dy*dy
    dt = F * (dx2 * dy2) / (nu * (dx2 + dy2))
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = t[1] - t[0]
    
    u   = np.zeros((Nx+1, Ny+1))      # solution array
    u_n = np.zeros((Nx+1, Ny+1))      # solution at t-dt
     
# Set initial condition
    for i in range(0,Nx+1):
        for j in range (0,Ny+1):
         u_n[i,j] = I(x[i],y[j])
         
#Run through Nt timesteps    
    for n in range(0, Nt):       
      u[1:-1, 1:-1] = u_n[1:-1, 1:-1] + nu * dt * (
          (u_n[2:, 1:-1] - 2*u_n[1:-1, 1:-1] + u_n[:-2, 1:-1])/dx2
          + (u_n[1:-1, 2:] - 2*u_n[1:-1, 1:-1] + u_n[1:-1, :-2])/dy2 )
      
# Switch variables before next step
    u_n = u.copy()
    
    t1 = time.clock()
 
    return t1-t0,u_n,u 



