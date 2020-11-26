# -*- coding: utf-8 -*-
"""
author: Silvia Vargas

"""
import scipy as sp
"""
def diffuse(I, nu, L, F, T): #tolto dt come parametro
    import time;  t0 = time.clock()  # For measuring the CPU time
    dx,dy= L* 0.01
    dx2, dy2 = dx*dx, dy*dy
    dt = F * (dx2 * dy2) / (nu * (dx2 + dy2))
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    Nx,Ny = 100
    x,y = np.linspace(0, L, 101)       # Mesh points in space         
    
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
"""
class Diffusion(object):
    """
    Class which implements a numerical solution of the 2d diffusion equation
    """
    def __init__(self, dx, dy, nu, kind, nt, L):
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
                 self.u,self.ui = self.get_initial_conditions(kind)
                 
    def get_initial_conditions(self, kind):
        # Start u and ui off as zero matrices:
        ui = sp.zeros([self.nx,self.ny])
        u = sp.zeros([self.nx,self.ny])
        # Now, set the initial conditions (ui).
        for i in range(self.nx):
            for j in range(self.ny):
                if kind == "circle":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= self.L*0.03 ):
                        ui[i,j] = 1
                if kind == "two_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    q = (i*self.dx - self.L*0.2)**2 + (j*self.dy - self.L*0.2)**2
                    if ( p  <= self.L*0.03 or q <= self.L*0.03 ):
                        ui[i,j] = 1 
                elif kind == "square":
                    if ( i>int(0.4*self.nx) and i<int(0.6*self.nx) and j>int(0.4*self.nx) and j<int(0.6*self.nx)):
                        ui[i,j] = 1
                elif kind == "donut":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= .03*self.L and p >= 0.020*self.L ):
                        ui[i,j] = 1
                elif kind == "concentric_circles":
                    p = (i*self.dx - self.L*0.5)**2 + (j*self.dy - self.L*0.5)**2
                    if ( p  <= .03*self.L and p >= 0.020*self.L ):
                        ui[i,j] = 1
                    elif ( p  <= .108*self.L and p >= 0.09*self.L ):
                        ui[i,j] = 1
                elif kind == "rod" :
                    if ( i>int(0.4*self.nx) and i<int(0.45*self.nx) and j>int(0.4*self.ny) and j<int(0.7*self.ny)):
                        ui[i,j] = 1
        return u,ui
    
    def evolve_ts(self):
        self.u[1:-1, 1:-1] = self.ui[1:-1, 1:-1] + self.nu*self.dt*( (self.ui[2:, 1:-1] - 2*self.ui[1:-1, 1:-1] + self.ui[:-2, 1:-1])/self.dx2 + (self.ui[1:-1, 2:] - 2*self.ui[1:-1, 1:-1] + self.ui[1:-1, :-2])/self.dy2 )
        self.ui = self.u.copy()
        

         
               
                



