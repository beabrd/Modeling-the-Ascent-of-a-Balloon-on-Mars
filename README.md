# Modeling-the-Ascent-of-a-Balloon-on-Mars

Welcome to our code!

All of the necessary libraries are imported at the start of each file. The easiest way to run the code is to paste it into a Google Colab notebook.

Each file produces plots for one numerical method. We coded Euler's method, improved Euler's method, and RK4 without the use of numerical integration libraries. The "rk4_5" file uses the solve_ivp function from the SciPy library, applying a Runge_Kutta method of order 5. We used this to check the correctness of our code for the other numerical methods, but we included it here for comparison's sake. The two files titled "mass_vs_altitude" and "radius_vs_altitude" produce graphs that compare the altitudes reached by balloons of different masses and radii. 

In each file, it is possible to choose the simulation time and the timestep (h). It is also possible to change the initial velocity assumption. These changes can be made by modifying the initial value of the variables in the section of the code titled "inital conditions: modifiable" in each file.

For short-duration simulations and zoomed-in analysis, smaller time steps improve accuracy and allow clearer comparison between numerical methods. For long-duration simulations (simulation times greater than 1000 s), a time step of h=1s or smaller provides the best balance between numerical stability, accuracy, and computational efficiency. 

Required libraries:
This code relies on standard scientific Python libraries, including NumPy for numerical computations and Matplotlib for data visualization. The SciPy library is also required for the rk4_5 file, which uses the solve_ivp adaptive Runge–Kutta solver. Please ensure that all required libraries are installed before running the code.

