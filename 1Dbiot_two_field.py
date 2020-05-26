#!/usr/bin/env python3

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def E_nu_to_mu_lmbda(E, nu):
    mu = E/(2*(1.0+nu))
    lmbda = (nu*E)/((1-2*nu)*(1+nu))
    return (mu, lmbda)

L = 15 # mm
N = 60
dt = 0.1
mesh = IntervalMesh(N, 0.0, L) 


#mesh = UnitSquareMesh(n, n)
cell = mesh.ufl_cell()
dim = mesh.topology().dim()

V = VectorElement("CG", cell, 1)
Q = FiniteElement("CG", cell, 1)

W = FunctionSpace(mesh, V*Q) 

x0 = W.tabulate_dof_coordinates()[:,0]
indices = np.argsort(x0)
x0 = x0[indices]

# Material parameters
E = 584.0 # Pa 
nu = 0.35 
(mu, lmbda) = E_nu_to_mu_lmbda(E, nu)

def sigma_star(u):
    I = Identity(mesh.topology().dim())
    return 2*mu*sym(grad(u)) + lmbda*div(u)*I


#(u, p) = TrialFunctions(W)
w = Function(W)
u, p = split(w)
(v, q) = TestFunctions(W)

# Solutions at previous timestep
w_ = Function(W)
u_, p_ = split(w_)

w_.vector()[:] = 0
# Material parameters
c = Constant(0.0)#10.0**(-10))
#c = Constant(0.000001)
alpha = Constant(1.0)
K = Constant(1.57*10**(-5))
# K = Constant(1.0*10**(-5))

# Timestepping parameters
k = Constant(dt)


dpdt = (p - p_)/k
dudt = (u - u_)/k

theta = Constant(1.0)
um = theta*u + (1 - theta)*u_
pm = theta*p + (1 - theta)*p_

# Body sources
f = Constant((5.0,)*dim)
g = Constant(0.0)

# top 
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 15.0)
    
# bottom
class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 15.0)    
        
top = top()
bottom = bottom()
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
top.mark(boundary_subdomains,1) 
bottom.mark(boundary_subdomains,2) 

# Essential boundary conditions
u_bottom = Constant((0.0,)*dim)
p_top = Constant(0.0)
#p_top = Expression("A*(1.0 - x[0]/L) + B", degree=1, A=13.0, B=0.0, L=L)
# Dirichlet BC
# Top -> zero pressure
bc0 = DirichletBC(W.sub(1), p_top, "near(x[0], 15.0) && on_boundary")
# Bottom -> zero displacement
bc1= DirichletBC(W.sub(0), u_bottom, "near(x[0], 0.0) && on_boundary")

bcs = [bc0, bc1]

# Neumann Bc
# Top(0) -> Load
# Bottom(1) -> no flow -> no need to do anything 

# Boundary sources
n = FacetNormal(mesh)
#t = - alpha*p_bar*n
t = as_vector((-1000.0,))
s = Constant(0.0)


dss = ds(subdomain_data =boundary_subdomains)
# Variational formulation
F = (inner(sigma_star(um), grad(v)) + alpha*pm*div(v) - inner(f, v))*dx - inner(t, v)*dss(1) \
     +(- c*dpdt*q + alpha*div(dudt)*q - inner(K*grad(pm), grad(q)) - g*q)*dx + s*q*dss(1)


#L0 = (inner(sigma_star(um), grad(v)) - alpha*pm*div(v) - inner(f, v))*dx - inner(t,v)*dss(1)
#L1 = (-c*dpdt*q - alpha*div(dudt)*q - inner(K*grad(pm), grad(q)) + g*q)*dx 

#F = L0 +L1


#(a, L) = system(F)



time = 0.0
count = 0
T = 9000.0
file_p = File("biot2p-nonlinear/p.pvd", "compressed")
file_u = File("biot2p-nonlinear/u.pvd", "compressed")

while (time < (T + DOLFIN_EPS)):
    
    solve(F == 0, w, bcs)
    if count % 1000 == 0:
        print("Time: ", time)
        (u_, p_) = w.split(deepcopy=True)
        #plt.plot(p_.vector(), label = str(time))
        file_p << p_
        #file_u << p_
    #plot(p, key="p", title="Pressure")
    #plot(-K*grad(p), key="v", title="Fluid velocity")
    
    # Update solutions
    time = time + dt
    w_.assign(w)
    count = count + 1
