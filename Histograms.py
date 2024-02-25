#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:11:20 2021

Plots

@author: liz
"""

sim = Simulation (100, 30, 0.5)
sim.run(50, True)

#print(sim._total_momentum)

# print(sim.temp())
# print(sim.pressure())
# print( sim._duration, sim._impulse, sim.pressure())

#%%

"""
Plot histograms of 
- distance to centre for all balls
- separation of balls

"""

dist_cent = [] #list of distance to centre for all balls
for i in range(1, len(sim._ball)-1):
    dist = np.sqrt(abs(np.dot(sim._ball[i+1]._pos,sim._ball[i+1]._pos)))
    dist_cent.append(dist)

f2 = pl.figure(figsize = [5,5])
y, binEdges = np.histogram(dist_cent)
plt.hist(dist_cent, color = "slateblue")
bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])

#fit linear plot to bin centres
fit1, cov1 = np.polyfit(bincenters, y,1, cov = True)
pfit1 = np.poly1d(fit1)

plt.plot(bincenters, y, 'o', color='black')
plt.plot(bincenters, pfit1(bincenters), color="black")
plt.grid()
plt.figtext(0,0,s="R = 30, r = 0.5, num = 100, collisions = 100")
plt.xlabel("Distance to centre")
plt.ylabel("Frequency")
plt.title("Histogram of distance of balls from centre")
plt.savefig("Histogram of distance of balls from centre")
plt.show()

#%%
ball_sep = [] #list of separations for all balls 
for i in range(1, len(sim._ball)-1):
    for j in range(1, len(sim._ball)-1):
        if i > j:
            sep = np.sqrt(abs(np.dot(sim._ball[i]._pos, sim._ball[j]._pos)))
            ball_sep.append(sep)
            
f3 = pl.figure(figsize = [5,5])
y1, binEdges1 = np.histogram(ball_sep, bins = 50)
bincenters1 = 0.5 * (binEdges1[1:] + binEdges1[:-1])

def ball_sep_func(r, R, A):
    return A*(4*r/(np.pi*R*R))*(np.arccos(r/(2*R))-(r/(2*R))*np.sqrt((1-(r*r)/(4*R*R))))
                         
    
params, cov = optimize.curve_fit(ball_sep_func, bincenters1, y1, p0=[15, 10e3])

plt.hist(ball_sep,bins=50, color = "darkslateblue")
plt.plot(bincenters1, y1, "o", color="black")
plt.plot(bincenters1, ball_sep_func(bincenters1, params[0], params[1]), "-",color = "black")
plt.grid()
plt.figtext(0,0,s="R = 30, r = 0.5, num = 100, collisions = 100")
plt.xlabel("Distance between balls")
plt.ylabel("Frequency")
plt.title("Histogram of distance between balls")
plt.savefig("Histogram of distance between balls")
plt.plot()

print("R = ", params[0], "A = ", params[1])

#%%

#Kinetic energy vs time plot

#print(sim._total_ke)

f4 = pl.figure(figsize = [5,5])
plt.plot(sim._intervals, sim._total_ke, "o", color = "mediumslateblue")
plt.grid()
plt.figtext(0,0,s="R = 30, r = 0.5, num = 100, collisions = 100")
plt.xlabel("Time, s")
plt.ylabel("Total kinetic energy of system, J")
plt.title("Kinetic energy of the system vs time")
plt.savefig("Kinetic energy of the system vs time")
plt.show()
#%%

#Momentum vs time plot

f5 = pl.figure(figsize = [5,5])
mom_mag = [] 
for i in range (len(sim._total_momentum)):
    mag = np.sqrt(abs(np.dot(sim._total_momentum[i], sim._total_momentum[i])))
    mom_mag.append(mag)
plt.plot(sim._intervals, mom_mag,"o", color="mediumpurple")
plt.grid()
plt.figtext(0,0,s="R = 30, r = 0.5, num = 100, collisions = 100")
plt.xlabel("Time, s")
plt.ylabel("Total momentum of the system")
plt.title("Momentum of the system vs time")
plt.savefig("Momentum energy of the system vs time")
plt.show()

#%%

#Pressure vs temperature plot

#this info was found manually by changing velocities and printing temp/pressure

temps = [8.227744217969322e+22, 2.2656567663908656e+22, 3.776002331074309e+22, 5.1877527431422685e+22, 7.339710056270283e+22, 1.1257922444317322e+23, 1.46826906400501e+23, 1.8762273476440028e+23, 2.135239351878287e+23, 2.588551759623009e+23]
pressures = [0.15724977094764642, 0.023336759292058393, 0.037693534733659086, 0.06793862956531824, 0.1431241723518843, 0.2402097827204003, 0.29759327831366483, 0.47678721872389485, 0.5621173646218258, 0.7108241495309895]

temps.sort()
pressures.sort()

#fit a linear graph on the data points
fit, cov = np.polyfit(temps, pressures,1, cov = True)
pfit = np.poly1d(fit)

plt.plot(temps, pressures,"o", color = "rebeccapurple")
plt.plot(temps, pfit(temps), color = "blueviolet")
plt.grid()
plt.figtext(0,0,s="R = 30, r = 0.5, num = 100, collisions = 100, data points = 10")
plt.xlabel("Temperature, K")
plt.ylabel("Pressure, Pa")
plt.title("Pressure vs temperature of the system")
plt.savefig("Pressure vs temperature of the system")
plt.show()

#gradient should be Nk/V
N = 100
k = 1.38e-23
V = np.pi * 30*30

theoretical_grad = N*k/V

gradient = fit[0]
gradient_unc = np.sqrt(cov[0,0])

print (gradient,  "+/-", gradient_unc)
print (theoretical_grad)

#%%
#Ball radius affecting PV=NkT
rads = [2, 1.5, 1, 0.5, 0.25, 0.125]

# temps1 = []
# pressures1 = []

# for i in range(len(rads)):
#     sim = Simulation (100, 30, rads[i])
#     sim.run(100, True)
#     temps1.append(sim.temp())
#     pressures1.append(sim.pressure())


# m = []
# for i in range (len(temps1)):
#     n = pressures1[i]/temps1[i]
#     m.append(n)

# theoretical_grad_array = []
# for i in range (len(m)):
#     theoretical_grad_array.append(theoretical_grad)

plt.plot (rads, m, "o", label = "Gradient of P-V plot", color="darkorchid")
plt.plot(rads, theoretical_grad_array, label = "Theoretical value of Nk/V", color = "indigo")
plt.grid()
plt.figtext(0,0,s="R = 30, num = 100, collisions = 100, data points = 6")
plt.xlabel("Radius of balls")
plt.ylabel("Value of gradient of the P-T plot")
plt.title("How changing the radius of the balls affects the equation of state")
plt.savefig("How changing the radius of the balls affects the equation of state")
plt.show()



#%%
"""
Maxwell-Boltzmann Distribution

"""

speeds = []
theoretical_mb = []

# Theoretical Maxwell-Boltzmann
# def mb (mass, temp, speed):
#     k = 1.38e-23
#     f = (mass/(np.pi*k*temp))**(3/2)*4*np.pi*speed*speed*math.exp((-mass/speed*speed)/2*k*temp)
#     return f

for i in range (len(sim._total_velocities)):
    speed = np.sqrt(abs(np.dot(sim._total_velocities[i], sim._total_velocities[i])))
    #mb_val = mb(sim._ball[1].m, sim.temp(), speed)
    speeds.append(speed)
    #theoretical_mb.append(mb_val)

y2, binEdges2 = np.histogram(speeds, bins = 10)
bincenters2 = 0.5 * (binEdges2[1:] + binEdges2[:-1])

def mb(v, m, T, A):
    return (v * A*np.exp(-0.5*m*v*v/(T*1.38e-23)))
                         
    
params2, cov2 = optimize.curve_fit(mb, bincenters2, y2, p0=[1, 2e23, 1500])

plt.hist(speeds, bins = 10, color = "thistle")
plt.plot(bincenters2, y2, "o", color="black")
plt.plot(bincenters2, mb(bincenters2, params2[0], params2[1], params2[2]), "-",color = "black")
plt.figtext(0,0,s="R = 30, r = 0.5 num = 100, collisions = 500")
plt.xlabel("Speed of balls")
plt.ylabel("Frequency")
plt.grid()
plt.title("Histogram of ball speeds")
plt.savefig("Histogram of ball speeds")

plt.show()

print (params2[0], params2[1], params2[2])


#%%
"""
Van der Waals
"""

def vdw (T, N, V, a, b):
    return ((N*1.38e-23)/(V-N*b))*T-a*(N/V)*(N/V)

temps = np.array([8.227744217969322e+22, 2.2656567663908656e+22, 3.776002331074309e+22, 5.1877527431422685e+22, 7.339710056270283e+22, 1.1257922444317322e+23, 1.46826906400501e+23, 1.8762273476440028e+23, 2.135239351878287e+23, 2.588551759623009e+23])
pressures = np.array([0.15724977094764642, 0.023336759292058393, 0.037693534733659086, 0.06793862956531824, 0.1431241723518843, 0.2402097827204003, 0.29759327831366483, 0.47678721872389485, 0.5621173646218258, 0.7108241495309895])

params3, cov3 = optimize.curve_fit(vdw, temps, pressures, p0=[100, 30*30*np.pi, 3e-6, 10])
plt.plot(temps, vdw(temps, params3[0], params3[1], params3[2], params3[3]), "-", color = "mediumvioletred", label = "Van der Waals fit")
plt.plot(temps, pressures, "o", color = "orchid")
plt.figtext(0,0,s="R = 30, r = 0.5 num = 100, collisions = 100")
plt.xlabel("Temperature, K")
plt.ylabel("Pressure, Pa")
plt.title("Plot of data used to fit to van der Waal’s law")
plt.grid()
plt.legend()
plt.savefig("Plot of data used to fit to van der Waal’s law")

plt.show()
print(params3[0], params3[1], params3[2], params3[3])
