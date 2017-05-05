import numpy as np
import matplotlib.pyplot as pt

mesh = np.loadtxt("OutputData/x_mesh.dat") # Holds the data of the mesh grid.
arrayh = np.loadtxt("OutputData/output.dat") # Holds the data of the computed solution
array = np.loadtxt("OutputData/true.dat"); # Holds the data of the true solution
h, err = np.loadtxt("OutputData/errors.dat", unpack=True)

ts = arrayh[:, 0] # time steps are in first column of output.dat file
us_h = arrayh[:,1:] # computed solutions are on columns: from 1 onwards, of file output.dat
us_true = array[:,1:] # true solutions are on columns: from 1 onwards, of file true.dat

_maph = {}
_map = {}

for i in range(ts.size):
    t = ts[i]
    _maph[t] = us_h[i] # Map time step to its computed solution.
    _map[t] = us_true[i] # Map time step to its true solution.

ts_toplot = [0, 1/2, 1] # Will only plot for this time values.

for time in ts_toplot:
    fig = pt.figure(figsize=(9,7))
    ax = fig.add_subplot(1,1,1)
    approx_sol = _maph[time] 
    true_sol = _map[time]
    ax.plot(mesh, approx_sol, 'ro', label = 'computed')
    ax.plot(mesh, true_sol, 'b-', label = 'true')
    ax.set_title("Case 2 \n Heat Equation @ time: "+str(time))
    ax.legend(loc='upper left')
    ax.set_xlabel('x space')
    ax.set_ylabel('u(x,'+str(time)+') state')
    fig.savefig("P2_time_"+str(time)+".png")

fig2 = pt.figure(figsize=(9,9))
ax2 = fig2.add_subplot(1,1,1)
ax2.loglog(h,err,'b-', lw=2, label = r'$errors$')
ax2.loglog(h, (h*h), 'ko', lw=2, label = r'$O(h^2)$')
ax2.legend(loc='upper left')
ax2.set_xlabel('Step Size')
ax2.set_ylabel('Error')
ax2.set_title("Case 2 \n "+r"Convergence of Implicit Method taking $\Delta t = 4h^2$"+"\nError values v.s. Step Size")
fig2.savefig("P2_errors.png")

pt.show()