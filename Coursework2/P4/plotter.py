import numpy as np
import matplotlib.pyplot as pt

mesh = np.loadtxt("OutputData/mesh.dat")
us = np.loadtxt("OutputData/computed.dat")
true_us = np.loadtxt("OutputData/unconstrained_true.dat")
array = np.loadtxt("OutputData/iterations.dat")

iters = array[:, 0] # iteration number is in first column of iterations.dat file
result = array[:,1:] # result are on columns: from 1 onwards, of file iterations.dat

_map = {}

for i in range(iters.size):
    _iter = iters[i]
    _map[_iter] = result[i] # Map iteration to its true solution.

# First 8 iterations
fig = pt.figure(figsize=(9,7))
ax = fig.add_subplot(1,1,1)
for k in range(1, 9):	
	approx_sol = _map[k] # result at iteration k
	to_plot = np.zeros(len(mesh))
	to_plot[1:-1] = approx_sol # BC's are zero at end points
	ax.plot(mesh, to_plot, label = "iteration"+str(k))
ax.set_title("Case 4 \n Iterations of the Projected SOR method $(\omega = 1.8)$")
ax.legend(loc='upper right', ncol=3)
ax.set_xlabel('x Space')
ax.set_ylabel('Approximate solution U(x)')
ax.set_ylim(0,1.5)
fig.savefig("P4_iterations.png")


fig1 = pt.figure(figsize=(9,7))
fig1.subplots_adjust(hspace=.5)
ax = fig1.add_subplot(1,1,1)
ax.plot(mesh, us, 'r-', label = 'Constrained Approximation')
ax.plot(mesh, true_us, 'b-', label = 'Unconstrained True')
ax.set_title("Case 4 \n Converged Approximation v.s. Unconstrained Approximation")
ax.legend(loc='lower center')
ax.set_xlabel('x space')
ax.set_ylabel('u(x) state')
fig1.savefig("P4_constrained_vs_unconstrained.png")

pt.show()