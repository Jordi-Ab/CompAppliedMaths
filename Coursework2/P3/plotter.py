import numpy as np
import matplotlib.pyplot as pt

case = 2

mesh = np.loadtxt("OutputData/grid_case"+str(case)+".dat") # Holds the data of the mesh grid.
arrayh = np.loadtxt("OutputData/computed_case"+str(case)+".dat") # Holds the data of the computed solution
array = np.loadtxt("OutputData/true_case"+str(case)+".dat"); # Holds the data of the true solution
h, err = np.loadtxt("OutputData/errors_case"+str(case)+".dat", unpack=True)

ts = arrayh[:, 0] # time steps are in first column of output.dat file
us_h = arrayh[:,1:] # computed solutions are on columns: from 1 onwards, of file output.dat
us_true = array[:,1:] # true solutions are on columns: from 1 onwards, of file true.dat

_maph = {}
_map = {}

for i in range(ts.size):
    t = ts[i]
    _maph[t] = us_h[i] # Map time step to its computed solution.
    _map[t] = us_true[i] # Map time step to its true solution.

ts_toplot = [0, 2.5, 5] # Will only plot for this time values.

for time in ts_toplot:
    fig = pt.figure(figsize=(9,7))
    ax = fig.add_subplot(1,1,1)
    approx_sol = _maph[time] 
    true_sol = _map[time]
    ax.plot(mesh, approx_sol, 'ro', label = 'computed')
    ax.plot(mesh, true_sol, 'b-', label = 'true')
    ax.set_title("Case 3("+"i"*case+") \n Black-Scholes @ time: "+str(time))
    ax.legend(loc='upper right')
    ax.set_xlabel('x \n Asset Price')
    ax.set_ylabel('Value of European Put Option \n v(x,'+str(time)+') ')
    fig.savefig("BS_Case"+str(case)+"_time_"+str(time)+".png")

fig2 = pt.figure(figsize=(9,9))
ax2 = fig2.add_subplot(1,1,1)
ax2.loglog(h,err,'b-', lw=2, label = r'$errors$')
ax2.loglog(h, (h**2), 'ko', lw=2, label = r'$O(h^2)$')
ax2.legend(loc='upper left')
ax2.set_xlabel('Step Size')
ax2.set_ylabel('Infinity Norm Error')
ax2.set_title(r"Case 3("+"i"*case+") \n Convergence of Semi Implicit Method taking $\Delta t = 4h^2$"+"\nError values v.s. Step Size")
fig2.savefig("BS_errors_case"+str(case)+".png")

pt.show()