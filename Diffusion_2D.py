# -*- coding: utf-8 -*-
"""
author: Silvia Vargas

"""
import scipy as sp
import matplotlib.pyplot as plt

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
                 self.nx = int(L/dx)
                 self.ny = int(L/dy)
                # For stability, this is the largest interval possible
                # for the size of the time-step:
                 self.dt = self.dx2*self.dy2/( 2*nu*(self.dx2+self.dy2) )
                # unknown u at new time level t and u at the previous 
                # time level t-dt
                 self.ui = self.get_initial_conditions(kind)
                 
    def get_initial_conditions(self, kind):#how to improve the implementation of IC/vectorization
        """Get the possible initial condition to solve the diffusion function, the options are: 
            "circle"
            "two_circles"
            "square"
            "donut"
            "concentric_circles"
            "rod"
            "semicircle"
        """
        # Start u and ui off as zero matrices:
        ui = sp.zeros([self.nx,self.ny])
        # Now, set the initial conditions (ui).
        for i in range(self.nx):
            for j in range(self.ny):
                if kind == "circle":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= self.L*0.1 ):
                        ui[i,j] = 1
                if kind == "two_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.25)**2
                    q = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.75)**2
                    if ( p  <= self.L*0.03 or q <= self.L*0.03 ):
                        ui[i,j] = 1 
                elif kind == "square":
                    if ( i>int(0.4*self.nx) and i<int(0.6*self.nx) and j>int(0.4*self.nx) and j<int(0.6*self.nx)):
                        ui[i,j] = 1
                elif kind == "donut":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= 0.5*self.L and p >= 0.3*self.L ):
                        ui[i,j] = 1
                elif kind == "concentric_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= .3*self.L and p >= 0.2*self.L ):
                        ui[i,j] = 1
                    elif ( p  <= .9*self.L and p >= 0.8*self.L ):
                        ui[i,j] = 1
                elif kind == "rod" :
                    if ( i>int(0.4*self.nx) and i<int(0.45*self.nx) and j>int(0.1*self.ny) and j<int(0.9*self.ny)):
                        ui[i,j] = 1
                elif kind == "semicircle" :
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= self.L*0.2 and j>=int(0.5*self.ny) ):
                        ui[i,j] = 1      
        return ui
    
    def evolve_ts(self):
        for i in range(self.nt + 1):#the cicle is nedeed? 
         self.u[1:-1, 1:-1] = self.ui[1:-1, 1:-1] + self.nu*self.dt*( (self.ui[2:, 1:-1]
                              - 2*self.ui[1:-1, 1:-1] + self.ui[:-2, 1:-1])/self.dx2 + 
                              (self.ui[1:-1, 2:] - 2*self.ui[1:-1, 1:-1] + self.ui[1:-1, :-2])/self.dy2 )
         self.ui = self.u.copy()
         
         return self.u, self.ui
        

         
#%%               
Diff1=Diffusion(0.1, 0.1, 1,"semicircle",100,10)     
ui=Diff1.get_initial_conditions("semicircle")          


plt.figure(figsize=(15, 7))


# plt.grid(True, which='major', linestyle='--', color='black', alpha=0.8)
plt.imshow(ui)
plt.show()


