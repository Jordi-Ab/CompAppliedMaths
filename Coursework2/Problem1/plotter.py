import numpy as np
import matplotlib.pyplot as pt

mesh = np.loadtxt("OutputData/mesh.dat")
us = np.loadtxt("OutputData/computed.dat")
true_us = np.loadtxt("OutputData/true.dat")
h, err = np.loadtxt("OutputData/errors.dat", unpack=True)

fig = pt.figure(figsize=(9,7))
fig.subplots_adjust(hspace=.5)
ax = fig.add_subplot(2,1,1)
ax.plot(mesh, us, 'ro', label = 'computed')
ax.plot(mesh, true_us, 'b-', label = 'true')
ax.set_title(r"Model Problem 1")
ax.legend(loc='upper left')
ax.set_xlabel('x space')
ax.set_ylabel('u(x) state')

ax2 = fig.add_subplot(2,1,2)
ax2.loglog(h,err,'b-', lw=2)
ax2.set_xlabel('Step Size')
#ax2.set_xlim(0, 1)
ax2.set_ylabel('Error')
ax2.set_title("Errors v.s. Step Size")

fig.savefig("Aprox_sol.png")
pt.show()