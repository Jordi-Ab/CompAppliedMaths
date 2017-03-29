import numpy as np
import matplotlib.pyplot as pt
import scipy.linalg as la

def grid_2norm(xvec, h):
    s = 0
    for x in xvec:
        s += abs(x)**2
    return np.sqrt(h*s)

def f(xvec): 
	return (np.pi**2)*np.sin(np.pi*xvec)

def true_sol(x):
    return x - np.sin(np.pi*x)

def solve(l_BC, r_BC, this_mesh):
    
    _n = this_mesh.size-1
    
    # Assume equally spaced points
    _h = this_mesh[1] - this_mesh[0]
    
    # Inner Second Order Difference Matrix
    first_row = np.zeros(_n-1)
    first_row[0] = -2; 
    first_row[1] = 1;
    A = (1/_h**2)*la.toeplitz(first_row)
    
    # Inner right hand side
    rhs = f(this_mesh[1:-1])
    
    # Apply Boundary Conditions to rhs
    rhs[0] = rhs[0] - l_BC/(_h**2)
    rhs[-1] = rhs[-1] - r_BC/(_h**2)
    
    # Solve inner system
    us = la.solve(A, rhs)
    
    # Assemble complete solution
    sol = np.zeros(_n+1)
    sol[0] = l_BC # Left Boundary condition
    sol[1:-1] = us # Computed result
    sol[-1] = r_BC # Right Boundary condition
    
    return sol

# MAIN

n = 200
alpha = 0
beta = 1

mesh = np.linspace(0, 1, n+1)
h = mesh[1] - mesh[0]

sol=solve(alpha, beta, mesh)

# Error
tr_sol = true_sol(mesh)
diff = tr_sol - sol
error = grid_2norm(diff, h)
print("For h: "+str(h)+",\t\t error="+str(error))

fig = pt.figure()
ax1 = fig.add_subplot(2,1,1)
ax1.plot(mesh, sol, 'ro', label='computed')
ax1.plot(mesh, true_sol(mesh), label='true')
ax1.set_title(r"Model Problem 1")
ax1.legend(loc='upper left')
ax1.set_xlabel('time')
ax1.set_ylabel('state')

# CODE VERIFICATION

print("Code Verification:")
mesh_size = 1/2
hs = []
errs = []
while(mesh_size>1e-3):
    mesh_size /= 2
    this_n = 1/mesh_size
    this_mesh = np.linspace(0,1,this_n)
    this_sol = solve(alpha, beta, this_mesh)
    true_s = true_sol(this_mesh)
    this_diff = true_s - this_sol
    #norm = grid_2norm(this_diff, mesh_size)
    norm = la.norm(this_diff, np.inf)
    
    hs.append(mesh_size)
    errs.append(norm)
    
    print("h = "+str(mesh_size)+ ", n = "+str(this_n)+"\t\t error="+str(norm))

ax2 = fig.add_subplot(2,1,2)
ax2.loglog(hs, errs, 'b-', lw=2)

pt.show()

