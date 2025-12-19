import numpy as np
import matplotlib.pyplot as plt


mass_total = 20  #total mass of ballon system (ballon + gas + payload) in kg
mass_gas = 5 #mass of lifting gas inside the ballon (kg)
g = 3.71 #gravitational acceleration on Mars in m/s^2

C_D = 0.47 #drag coefficient


R_gas = 2077.0 #specfic gas constant for gas **insert gas tyep** (**insert units)

r_max = 10
V_max = (4/3) * np.pi *r_max**3
#inital conditions
t0=0.00
tf = 10000
h=1

z0=0.0
v0=0.0

y0=np.array([z0,v0])

def mars_temperature_C(z):
    #return Martian temperature depending on alltiude
    if z > 7000: #z is altitude in meters
        return -23.4 - 0.00222 *z
    else:
        return -31.0 - 0.000998 * z

def mars_pressure(z):
    #return mars pressure in kpa depending on altitude
    return 0.699 * np.exp(-0.00009*z)

def mars_density (z):
    #return the density of martian air depending on altiude (z), tempeture, and pressure
    T_C = mars_temperature_C(z) #find tempeture in celesis
    T_K = T_C +273.15 #convert tempeture to kelvin
    p = mars_pressure(z) *1000 #find pressure and convert to pa
    return p/ (192.1 * T_K)

def balloon_volume(z):
    T_C = mars_temperature_C(z) #find tempeture in celesis
    T_K = T_C +273.15 #convert tempeture to kelvin
    p = mars_pressure(z)*1000  #find pressure
    V_ideal= (mass_gas*R_gas*T_K)/p
    return min(V_ideal, V_max)

def balloon_cross_sectional_area(z):
    V=balloon_volume(z) #cross sectional area for drag of sphere is pir^2
    return np.pi * ((3*V) / (4*np.pi)) **(2/3)

def derivatives(t,y):
    #returns derivatives of state variables y[0]=z(altitude) and y[1]= v(vertical velocity)
    #dy0/dt = v
    #dy1/dt = a
    #dy/dt = [v,a]

    z1=y[0]
    z = max(z1, 0.0) #insures balloon cant go underground
    v=y[1]

    rho=mars_density(z)
    V= balloon_volume(z)
    S = balloon_cross_sectional_area(z)

    #calulate buoyancy - weight term
    F_buoy = rho * V * g
    F_weight = mass_total * g


    #calulate drag term
    F_drag= 0.5*rho*C_D*S*v*abs(v) # allows drag to be negative

    #find net acceleration
    a = ((F_buoy - F_weight)-F_drag)/mass_total

    return np.array([v,a])


def RK4_solver(t0,tf,h,y0):
    #t0 is start time
    #tf is final time
    #h is time step
    #y0 is inital state[z0,v0]

    #np.arrange(start,stop,step)
    t_values = np.arange(t0,tf+h,h)

    # store solutions where rows time steps and colums is [z,v]
    y_values = np.zeros((len(t_values), len(y0)))

    #set the inital conditions, first column is y0
    y_values[0,:] = y0

    for n in range(len(t_values)-1):
        t_n = t_values[n]
        y_n = y_values[n,:]

        #compute slope f(t_n, y_n)
        k1 = derivatives(t_n,y_n)
        k2 = derivatives(t_n+h/2,y_n+(h/2)*k1)
        k3 = derivatives(t_n+h/2,y_n+(h/2)*k2)
        k4 = derivatives(t_n+h,y_n+h*k3)

        y_values[n+1,:] = y_n + (h/6)*(k1 +2*k2 +2*k3 +k4)


    return t_values, y_values




t1, y1 = RK4_solver(t0,tf,h,y0)

z1 = np.maximum(y1[:, 0], 0.0)
v1=y1[:,1]

plt.figure()
plt.plot(t1,z1)
plt.xlabel('Time (seconds)')
plt.ylabel('Altitude z (meters)')
plt.title('RK4_Solver: Altitude vs Time')
plt.grid(True) #do we want gridlines?
plt.show()

plt.figure()
plt.plot(t1,v1)
plt.xlabel('Time (seconds)')
plt.ylabel('Vertical Velocity v (m/s)')
plt.title('RK4_Solver: Velocity vs Time')
plt.grid(True) #do we want gridlines?
plt.show()
