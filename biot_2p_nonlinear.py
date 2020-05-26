#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 20:16:01 2019

@author: ranjeet
"""

from __future__ import division

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def E_nu_to_mu_lmbda(E, nu):
    mu = E/(2*(1.0+nu))
    lmbda = (nu*E)/((1-2*nu)*(1+nu))
    return (mu, lmbda)

# Create mesh
nx = 10
ny = 60
L = 15.0 # mm
N = 100

mesh =  RectangleMesh.create([Point(0.0,0.0), Point(1.0, L)],[nx,ny], CellType.Type.triangle)
n = FacetNormal(mesh)

#mesh = UnitSquareMesh(n, n)
cell = mesh.ufl_cell()
dim = mesh.topology().dim()

def sigma_star(u):
    I = Identity(mesh.topology().dim())
    return 2*mu*sym(grad(u)) + lmbda*div(u)*I

# top 
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 15.0)
    
# bottom
class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0)    

        
top = top()
bottom = bottom()
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
top.mark(boundary_subdomains,1) 
bottom.mark(boundary_subdomains,2) 

dt = 0.1

V = VectorElement("CG", cell, 2, 2)
Q = FiniteElement("CG", cell, 1)

W = FunctionSpace(mesh, V*Q) 

x0 = W.tabulate_dof_coordinates()[:,0]
indices = np.argsort(x0)
x0 = x0[indices]

# Material parameters
E = 584.0 # Pa 
nu = 0.35 
(mu, lmbda) = E_nu_to_mu_lmbda(E, nu)

(mu, lmbda) = 40e6, 40e6

# Material parameters
c = Constant(0.5)#10.0**(-10))
#c = Constant(0.000001)
alpha = Constant(1.0)
K = Constant(1.57*10**(-5))
# K = Constant(1.0*10**(-5))




#(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

# Solutions at previous timestep
w_ = Function(W)
u_, p_ = split(w_)

w = Function(W)
(u, p) = split(w)



# Body sources
f = Constant((0.0, 0.0))
g = Constant(0.0)
s = Constant(5.0)


# Dirichlet boundary conditions
u_bottom = Constant(0.0)
p_top = Constant(0.0)

# Bottom -> zero displacement
bc0= DirichletBC(W.sub(0).sub(1), u_bottom, "near(x[1], 0.0) && on_boundary") # No normal displacement for solid on bottom 
bc1 = DirichletBC(W.sub(0).sub(0), u_bottom, "on_boundary && (near(x[0], 0.0) || near(x[0], 1.0))") # No normal displacement for solid on sides
# Top -> zero pressure
bcp = DirichletBC(W.sub(1), p_top, "near(x[1], 15.0) && on_boundary")
bcs = [bc0, bc1, bcp]


# Neumann Bc
# Top(0) -> Load
t = Constant((0.0, -5))

dss = ds(subdomain_data =boundary_subdomains) 
# Timestepping parameters
k = Constant(dt)
dpdt = (p - p_)/k
dudt = (u - u_)/k

theta = Constant(0.5)
um = theta*u + (1 - theta)*u_
pm = theta*p + (1 - theta)*p_

## Variational formulation
#F = (inner(sigma_star(um), grad(v)) - alpha*pm*div(v) - inner(f, v)
#     - c*dpdt*q - alpha*div(dudt)*q - inner(K*grad(pm), grad(q)) + g*q)*dx \
#     + (inner(t, v) + s*q)*ds
##(a, L) = system(F)
#     
#     
## Variational formulation
#F = (inner(sigma_star(um), grad(v)) - alpha*pm*div(v) - inner(f, v)
#     - c*dpdt*q - alpha*div(dudt)*q - inner(K*grad(pm), grad(q)) + g*q)*dx \
#     + (inner(t, v) + s*q)*ds
     
# Variational formulation
F = (inner(sigma_star(um), grad(v)) + alpha*pm*div(v) - inner(f, v)
     - c*dpdt*q + alpha*div(dudt)*q - inner(K*grad(pm), grad(q)) - g*q)*dx \
     - (dot(t, v) - s*q)*ds
     

time = 0.0
count = 0
T = 50.0
file_p = File("biot2p-nonlinear/p.pvd", "compressed")
file_u = File("biot2p-nonlinear/u.pvd", "compressed")

while (time < (T + DOLFIN_EPS)):
    
    if count % 50 == 0:
        print("Time: ", time)
        (u1, p1) = w.split(deepcopy=True)
        #plt.plot(p_.vector(), label = str(time))
        file_p << p1
        file_u << u1
    
    solve(F == 0, w, bcs)
    
    # Update solutions
    time = time + dt
    w_.assign(w)
    count = count + 1
