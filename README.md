# The two-dimensional diffusion equation

The two-dimensional diffusion equation is:


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;\nu\frac{\partial^2u}{\partial&space;x^2}&plus;\nu\frac{\partial^2u}{\partial&space;y^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;\nu\frac{\partial^2u}{\partial&space;x^2}&plus;\nu\frac{\partial^2u}{\partial&space;y^2}" title="\frac{\partial u}{\partial t} = \nu\frac{\partial^2u}{\partial x^2}+\nu\frac{\partial^2u}{\partial y^2}" /></a>

where u is the quantity to calculate, *t* is for temporal variable, *x* and *y* are for spatial variables, and <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a> is the diffusivity constant.  

Before finding the solution *u* everywhere in *x* and *y*, and over time *t* we have to discretize the second derivative.

## Discretization of the second derivative

The second-order derivative can be represented geometrically as the line tangent to the curve given by the first derivative. We will discretize the second-order derivative with the Finite-difference method, it is a numerical method for solving differential equations by approximating derivative with finite differences: a combination of Forward Difference and Backward Difference of the first derivative. Consider the Taylor expansion of <a href="https://www.codecogs.com/eqnedit.php?latex=u_{i&plus;1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i&plus;1}" title="u_{i+1}" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=u_{i-1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i-1}" title="u_{i-1}" /></a> around <a href="https://www.codecogs.com/eqnedit.php?latex=u_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_i" title="u_i" /></a>:

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{i&plus;1}&space;=&space;u_i&space;&plus;&space;\Delta&space;x&space;\frac{\partial&space;u}{\partial&space;x}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^2}{2}&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^3}{3!}&space;\frac{\partial&space;^3&space;u}{\partial&space;x^3}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i&plus;1}&space;=&space;u_i&space;&plus;&space;\Delta&space;x&space;\frac{\partial&space;u}{\partial&space;x}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^2}{2}&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^3}{3!}&space;\frac{\partial&space;^3&space;u}{\partial&space;x^3}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" title="u_{i+1} = u_i + \Delta x \frac{\partial u}{\partial x}\bigg|_i + \frac{\Delta x^2}{2} \frac{\partial ^2 u}{\partial x^2}\bigg|_i + \frac{\Delta x^3}{3!} \frac{\partial ^3 u}{\partial x^3}\bigg|_i + O(\Delta x^4)" /></a>


<a href="https://www.codecogs.com/eqnedit.php?latex=u_{i-1}&space;=&space;u_i&space;-&space;\Delta&space;x&space;\frac{\partial&space;u}{\partial&space;x}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^2}{2}&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;-&space;\frac{\Delta&space;x^3}{3!}&space;\frac{\partial&space;^3&space;u}{\partial&space;x^3}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i-1}&space;=&space;u_i&space;-&space;\Delta&space;x&space;\frac{\partial&space;u}{\partial&space;x}\bigg|_i&space;&plus;&space;\frac{\Delta&space;x^2}{2}&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;-&space;\frac{\Delta&space;x^3}{3!}&space;\frac{\partial&space;^3&space;u}{\partial&space;x^3}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" title="u_{i-1} = u_i - \Delta x \frac{\partial u}{\partial x}\bigg|_i + \frac{\Delta x^2}{2} \frac{\partial ^2 u}{\partial x^2}\bigg|_i - \frac{\Delta x^3}{3!} \frac{\partial ^3 u}{\partial x^3}\bigg|_i + O(\Delta x^4)" /></a>


If we add these two expansions, you can see that the odd-numbered derivative terms will cancel each other out. If we neglect any terms of <a href="https://www.codecogs.com/eqnedit.php?latex=O(\Delta&space;x^4)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?O(\Delta&space;x^4)" title="O(\Delta x^4)" /></a> or higher (and really, those are very small), then we can rearrange the sum of these two expansions to solve for our second-derivative.

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{i&plus;1}&space;&plus;&space;u_{i-1}&space;=&space;2u_i&plus;\Delta&space;x^2&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i&plus;1}&space;&plus;&space;u_{i-1}&space;=&space;2u_i&plus;\Delta&space;x^2&space;\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i&space;&plus;&space;O(\Delta&space;x^4)" title="u_{i+1} + u_{i-1} = 2u_i+\Delta x^2 \frac{\partial ^2 u}{\partial x^2}\bigg|_i + O(\Delta x^4)" /></a>


Then rearrange to solve for <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^2&space;u}{\partial&space;x^2}\bigg|_i" title="\frac{\partial ^2 u}{\partial x^2}\bigg|_i" /></a> and the result is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;^2&space;u}{\partial&space;x^2}=\frac{u_{i&plus;1}-2u_{i}&plus;u_{i-1}}{\Delta&space;x^2}&space;&plus;&space;O(\Delta&space;x^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;^2&space;u}{\partial&space;x^2}=\frac{u_{i&plus;1}-2u_{i}&plus;u_{i-1}}{\Delta&space;x^2}&space;&plus;&space;O(\Delta&space;x^2)" title="\frac{\partial ^2 u}{\partial x^2}=\frac{u_{i+1}-2u_{i}+u_{i-1}}{\Delta x^2} + O(\Delta x^2)" /></a>

## Solution of the discretized 2D diffusion equation

We can now write the discretized version of the diffusion equation in 2D:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{u_{i,j}^{n&plus;1}&space;-&space;u_{i,j}^n}{\Delta&space;t}&space;=&space;\nu&space;\frac{u_{i&plus;1,j}^n&space;-&space;2&space;u_{i,j}^n&space;&plus;&space;u_{i-1,j}^n}{\Delta&space;x^2}&space;&plus;&space;\nu&space;\frac{u_{i,j&plus;1}^n-2&space;u_{i,j}^n&space;&plus;&space;u_{i,j-1}^n}{\Delta&space;y^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{u_{i,j}^{n&plus;1}&space;-&space;u_{i,j}^n}{\Delta&space;t}&space;=&space;\nu&space;\frac{u_{i&plus;1,j}^n&space;-&space;2&space;u_{i,j}^n&space;&plus;&space;u_{i-1,j}^n}{\Delta&space;x^2}&space;&plus;&space;\nu&space;\frac{u_{i,j&plus;1}^n-2&space;u_{i,j}^n&space;&plus;&space;u_{i,j-1}^n}{\Delta&space;y^2}" title="\frac{u_{i,j}^{n+1} - u_{i,j}^n}{\Delta t} = \nu \frac{u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n}{\Delta x^2} + \nu \frac{u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n}{\Delta y^2}" /></a>

We reorganize the discretized equation and solve for <a href="https://www.codecogs.com/eqnedit.php?latex=u_{i,j}^{n&plus;1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i,j}^{n&plus;1}" title="u_{i,j}^{n+1}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{i,j}^{n&plus;1}&space;=&space;u_{i,j}^n&space;&plus;&space;\frac{\nu&space;\Delta&space;t}{\Delta&space;x^2}(u_{i&plus;1,j}^n&space;-&space;2&space;u_{i,j}^n&space;&plus;&space;u_{i-1,j}^n)&space;&plus;&space;\frac{\nu&space;\Delta&space;t}{\Delta&space;y^2}(u_{i,j&plus;1}^n-2&space;u_{i,j}^n&space;&plus;&space;u_{i,j-1}^n)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{i,j}^{n&plus;1}&space;=&space;u_{i,j}^n&space;&plus;&space;\frac{\nu&space;\Delta&space;t}{\Delta&space;x^2}(u_{i&plus;1,j}^n&space;-&space;2&space;u_{i,j}^n&space;&plus;&space;u_{i-1,j}^n)&space;&plus;&space;\frac{\nu&space;\Delta&space;t}{\Delta&space;y^2}(u_{i,j&plus;1}^n-2&space;u_{i,j}^n&space;&plus;&space;u_{i,j-1}^n)" title="u_{i,j}^{n+1} = u_{i,j}^n + \frac{\nu \Delta t}{\Delta x^2}(u_{i+1,j}^n - 2 u_{i,j}^n + u_{i-1,j}^n) + \frac{\nu \Delta t}{\Delta y^2}(u_{i,j+1}^n-2 u_{i,j}^n + u_{i,j-1}^n)" /></a>

Besides the second derivative we discretized also the spatial domain and the time interval *x*, *y*, and *t*.

<a href="https://www.codecogs.com/eqnedit.php?latex=x_i=i\Delta&space;x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i=i\Delta&space;x" title="x_i=i\Delta x" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=y_j=&space;j\Delta&space;y" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_j=&space;j\Delta&space;y" title="y_j= j\Delta y" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=t_n=&space;n\Delta&space;t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t_n=&space;n\Delta&space;t" title="t_n= n\Delta t" /></a>

Where *i*, *j*, and *n* are the steps for each difference for *x*, *y*, and *t* respectively.
It can be shown that the maximum time step, <a href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;t" title="\Delta t" /></a> that we can allow without the process becoming unstable is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;t&space;=&space;\frac{1}{2\nu}\frac{(\Delta&space;x\Delta&space;y)^2}{(\Delta&space;x)^2&space;&plus;&space;(\Delta&space;y)^2}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;t&space;=&space;\frac{1}{2\nu}\frac{(\Delta&space;x\Delta&space;y)^2}{(\Delta&space;x)^2&space;&plus;&space;(\Delta&space;y)^2}." title="\Delta t = \frac{1}{2\nu}\frac{(\Delta x\Delta y)^2}{(\Delta x)^2 + (\Delta y)^2}." /></a>

## Structure of the project


The project is divided in two main blocks:

- The module [Diffusion_2D](https://github.com/silviavargas/Diffusion_2D/blob/master/Diffusion_2D.py) contains the **Class Diffusion** that and implements a numerical solution of the 2D diffusion equation, the parameters of the class are: the size interval to build the square mesh grid for the spacial coordinates *dx*, the diffusivity constant <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a>, the initial conditions *kind*, the number of timestep *nt* and the lenght of the mesh grid *L*. In this Class is also implemented the method **get_initial_conditions** that creates different initial conditions created with *for cycles* and returns the 2D array corresponding to the solution for the first timestep, the different kind are "circle", "two_circles", "square", "donut", "concentric_circles", "rod", "semicircle", "grid". Lastly is implemented the method **evolve_ts**, the algorithm evolve the system every timestep and returns a 3D array that contains the 2D diffusion equation solution for each timestep.

- The module [Diff_animation_3D](https://github.com/silviavargas/Diffusion_2D/blob/master/Diff_animation_3D.py) has been created to realize an animation using the matplotlib library, that show the process of diffusion over time, in this module there is the function **update_plot**, that takes the parameters *frame_number*, *zarray* and *plot* to create the 3D animation that is saved in GIF format.

Finally in the folder [GIFs](https://github.com/silviavargas/Diffusion_2D/tree/master/GIFs) are stored all the GIFs associated with the Configurations and in the folder [Test](https://github.com/silviavargas/Diffusion_2D/blob/master/Test.py) are present some oracle test to make sure the code is working as it should.   

These are the steps in order to start the program and to plot the results:
1) First, the user has to choose between the different configurations for the diffusion, contained in the folder [Configurations](https://github.com/silviavargas/Diffusion_2D/tree/master/Configurations) and eventually write a new one, using the same syntax of the other files in the folder; if the user wants to do so,
he has to specify the diffusion parameters (*dx*, <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a>, *kind*, *nt* and *L*) and also the name under which the GIF is saved. 
2) Then, to start the Diffusion the user has to launch the file [Diff_animation_3D](https://github.com/silviavargas/Diffusion_2D/blob/master/Diff_animation_3D.py) which imports its parameters from one of the file in [Configurations](https://github.com/silviavargas/Diffusion_2D/tree/master/Configurations) using ConfigParser library; there could be different types of configurations for the diffusion, depending on the size of the space interval, diffusivity constant, Initial conditiions, timestep or the name of the resulting animation, so the user has to specify the configuration he wants when launching the simulation file from the command line with the syntax ***"python Diff_animation_3D.py Configurations\name_of_the_configuration"*** (for example, circle.txt).

## Examples

To provide some examples, this is how the simulation of a given configuration looks like, for the configuration circle.txt, grid.txt and rod.txt:
![Diffusion_circle](https://github.com/silviavargas/Diffusion_2D/blob/master/GIFs/Diffusion_circle.gif)
![Diffusion_grid](https://github.com/silviavargas/Diffusion_2D/blob/master/GIFs/Diffusion_grid.gif)
![Diffusion_rod](https://github.com/silviavargas/Diffusion_2D/blob/master/GIFs/Diffusion_rod.gif)



 
    















