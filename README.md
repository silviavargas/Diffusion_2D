# The two-dimensional diffusion equation

The two-dimensional diffusion equation is:


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;\nu\frac{\partial^2u}{\partial&space;x^2}&plus;\nu\frac{\partial^2u}{\partial&space;x^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;\nu\frac{\partial^2u}{\partial&space;x^2}&plus;\nu\frac{\partial^2u}{\partial&space;x^2}" title="\frac{\partial u}{\partial t} = \nu\frac{\partial^2u}{\partial x^2}+\nu\frac{\partial^2u}{\partial x^2}" /></a>

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













