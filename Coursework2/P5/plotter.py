import numpy as np
import matplotlib.pyplot as pt

case = 1

mesh = np.loadtxt("OutputData/grid_case"+str(case)+".dat") # Holds the data of the mesh grid.
arrayh = np.loadtxt("OutputData/computed_case"+str(case)+".dat") # Holds the data of the computed solution
array2 = np.loadtxt("OutputData/computed_case"+str(case)+"_2.dat");

ts = arrayh[:, 0] # time steps are in first column of output.dat file
us_h = arrayh[:,1:] # computed solutions are on columns: from 1 onwards, of file output.dat
us2 = array2[:,1:]

_maph = {}
_maph2 = {}

for i in range(ts.size):
    t = ts[i]
    _maph[t] = us_h[i] # Map time step to its computed solution.
    _maph2[t] = us2[i]

ts_toplot = [0, 2.5, 5] # Will only plot for this time values.

for time in ts_toplot:
    fig = pt.figure(figsize=(9,7))
    ax = fig.add_subplot(1,1,1)
    pay_off = _maph2[0]
    approx_sol = _maph[time]
    approx_sol2 = _maph2[time] 
    ax.plot(mesh, approx_sol, 'r-', label = 'American Put')
    ax.plot(mesh, approx_sol2, 'b-', label = 'European Put')
    ax.plot(mesh, pay_off, 'g-', label = 'Pay off')
    ax.set_title("Case "+str(case+4)+"\n Black-Scholes @ time: "+str(time))
    ax.legend(loc='upper right')
    ax.set_xlabel('x \n Asset Price')
    ax.set_ylabel('Value of Put Option \n v(x,'+str(time)+') ')
    fig.savefig("BS_Case"+str(case)+"_time_"+str(time)+".png")

pt.show()