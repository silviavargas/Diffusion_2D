# -*- coding: utf-8 -*-
"""
Created on Wed Nov ‎11 ‏‎11:23:29 2020

@author: Silvia Vargas
"""
import numpy as np
#import matplotlib.pyplot as plt

class Diffusion(object):
    """
    Class which implements a numerical solution of the 2d diffusion equation
    """
    def __init__(self, dx, dy, nu, kind, nt, L):#how does the class works, the configuration is nedeed? 
                 self.L = L   
                 self.dx = dx # Interval size in x-direction.
                 self.dy = dy # Interval size in y-direction.
                 self.nu = nu # Diffusion constant.
                 self.nt = nt  #Number of time-steps to evolve system.
                 self.dx2 = dx**2
                 self.dy2 = dy**2
                 self.nx = int(round(L/dx))
                 self.ny = int(round(L/dy))
                # For stability, this is the largest interval possible
                # for the size of the time-step:
                 self.dt = self.dx2*self.dy2/( 2*nu*(self.dx2+self.dy2) )
                 self.u = self.get_initial_conditions(kind)
                 
    def get_initial_conditions(self, kind):#how to improve the implementation of IC/vectorization
        """Get the possible initial condition to solve the diffusion function, the options are: 
            "circle"
            "two_circles"
            "square"
            "donut"
            "concentric_circles"
            "rod"
            "semicircle"
            "grid"
        """
        # Start u:
        u = np.zeros([self.nx,self.ny,self.nt])
        # Now, set the initial conditions (ui).
        for i in range(self.nx):
            for j in range(self.ny):
                if kind == "circle":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= self.L*0.1 ):
                        u[i,j,0] = 1
                if kind == "two_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.25)**2
                    q = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.75)**2
                    if ( p  <= self.L*0.03 or q <= self.L*0.03 ):
                        u[i,j,0] = 1 
                elif kind == "square":
                    if ( i>int(round(0.4*self.nx)) and i<int(round(0.6*self.nx)) and j>int(round(0.4*self.nx)) and j<int(round(0.6*self.nx))):
                        u[i,j,0] = 1
                elif kind == "donut":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= 0.5*self.L and p >= 0.3*self.L ):
                        u[i,j,0] = 1
                elif kind == "concentric_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= .3*self.L and p >= 0.2*self.L ):
                        u[i,j,0] = 1
                    elif ( p  <= .9*self.L and p >= 0.8*self.L ):
                        u[i,j,0] = 1
                elif kind == "rod" :
                    if ( i>int(round(0.4*self.nx)) and i<int(round(0.45*self.nx)) and j>int(round(0.1*self.ny)) and j<int(round(0.9*self.ny))):
                        u[i,j,0] = 1
                elif kind == "semicircle" :
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= self.L*0.2 and j>=int(0.5*self.ny) ):
                        u[i,j,0] = 1
                elif kind == "grid":
                    if ( i>int(round(0.2*self.nx)) and i<int(round(0.3*self.nx)) or i>int(round(0.7*self.nx)) and i<int(round(0.8*self.nx))  or 
                        j>int(round(0.2*self.ny)) and j<int(round(0.3*self.ny)) or j>int(round(0.7*self.ny)) and j<int(round(0.8*self.ny))):
                        u[i,j,0] = 1       

        return u
    
    def evolve_ts(self):
     #   u = np.zeros([self.nx,self.ny,self.nt])
        for n in range(1, self.nt):       
            self.u[1:-1, 1:-1, n] = self.u[1:-1, 1:-1, n-1] + self.nu * self.dt * (
            (self.u[2:, 1:-1, n-1] - 2*self.u[1:-1, 1:-1, n-1] + self.u[:-2, 1:-1,n-1])/self.dx2
            + (self.u[1:-1, 2: ,n-1] - 2*self.u[1:-1, 1:-1, n-1] + self.u[1:-1, :-2, n-1])/self.dy2 )
            
         
        return self.u #, self.ui
        
"""
#%%               
Diff1=Diffusion(0.1, 0.1, 1,"grid",1000,10)     
#ui=Diff1.get_initial_conditions("two_circles")    #da togliere
Diff1.evolve_ts()
T=np.zeros([Diff1.nx, Diff1.ny])      
for i in range (0,Diff1.nx) :
    for j in range (0,Diff1.ny): 
        T[i,j]=Diff1.u[i,j,Diff1.nt-1]
     

plt.figure(figsize=(15, 7))


# plt.grid(True, which='major', linestyle='--', color='black', alpha=0.8)
plt.imshow(T)
plt.show()
"""
