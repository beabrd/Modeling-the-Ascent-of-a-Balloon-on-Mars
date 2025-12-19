from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


mass_total = 20  #total mass of ballon system (ballon + gas + payload) in kg
mass_gas = 5 #mass of lifting gas inside the ballon (kg)
g = 9 #gravitational acceleration on Mars in m/s^2

C_D = 0.47 #drag coefficient


R_gas = 2077.0 #specfic gas constant for Helium

r_max = 10
V_max = (4/3) * np.pi *r_max**3

t0=0
z0=0.0

#inital conditions: modifiable
tf = 10000
h=1
v0=0.0
###

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
    F_drag= -0.5*rho*C_D*S*v*abs(v) # allows drag to be negative

    #find net acceleration
    a = ((F_buoy - F_weight)+F_drag)/mass_total

    return np.array([v,a])


sol = solve_ivp(
    derivatives,
    t_span=(t0, tf),
    y0=y0,
    method='RK45',     # Runge–Kutta (4,5)
    max_step=h,        # force comparable timestep
    dense_output=False
)




t = sol.t
z = sol.y[0]
v = sol.y[1]

plt.figure()
plt.plot(t, z)
plt.xlabel('Time (s)')
plt.ylabel('Altitude z (m)')
plt.title('SciPy RK45: Altitude vs Time')
plt.grid(True)
plt.show()

plt.figure()
plt.plot(t, v)
plt.xlabel('Time (s)')
plt.ylabel('Velocity v (m/s)')
plt.title('SciPy RK45: Velocity vs Time')
plt.grid(True)
plt.show()
