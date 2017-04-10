import numpy as np
import matplotlib.pyplot as pt

case = 1

mesh = np.loadtxt("OutputData/mesh.dat")
us = np.loadtxt("OutputData/computed.dat")
true_us = np.loadtxt("OutputData/true.dat")
#h, err = np.loadtxt("OutputData/errors_case"+str(case)+".dat", unpack=True)

fig1 = pt.figure(figsize=(9,7))
fig1.subplots_adjust(hspace=.5)
ax = fig1.add_subplot(1,1,1)
ax.plot(mesh, us, 'ro', label = 'computed')
ax.plot(mesh, true_us, 'b-', label = 'true')
ax.set_title(r"Model Problem 1")
ax.legend(loc='upper left')
ax.set_xlabel('x space')
ax.set_ylabel('u(x) state')

#fig2 = pt.figure(figsize=(9,7))
#fig2.subplots_adjust(hspace=.5)
#ax2 = fig2.add_subplot(1,1,1)
#ax2.loglog(h,err,'b-', lw=2, label = r'$errors$')
#ax2.loglog(h, (h*h), 'ko', lw=2, label = r'$O(h^2)$')
#ax2.legend(loc='upper left')
#ax2.set_xlabel('Step Size')
#ax2.set_ylabel('Error')
#ax2.set_title("Errors v.s. Step Size for Case "+str(case))

fig1.savefig("P4_Uh_vs_U.png")
#fig2.savefig("P1_errors_case"+str(case)+".png")
pt.show()