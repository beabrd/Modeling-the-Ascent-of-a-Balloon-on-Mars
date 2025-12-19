# Modeling-the-Ascent-of-a-Balloon-on-Mars

Welcome to our code!

All of the necessary liraries are imported at the start of each file. 

Each file produces plots for one numerical method. We coded Euler's method, improved Euler's method, and RK4 without the use of numerical integration libraries. The "rk4_5" file uses the solve_ivp function from the SciPy library, applying a Runge_Kutta method of order 5. We used this to check the correctness of our code for the other numerical methods, but we included it here for comparison's sake. The two files titled "mass_vs_altitude" and "radius_vs_altitude" produce graphs that compare the altitudes reached by balloons of different masses and radii. 

In each file, it is possible to choose the simulation time and the timestep (h). It is also possile to change the initial velocity assumption. These changes can be made by modifying the initial value of the variables in the section of the code titled "inital conditions: modifiable" in each file.

IMPORTANT NOTES:
