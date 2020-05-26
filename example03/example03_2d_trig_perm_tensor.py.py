#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:43:08 2020

@author: ranjeet
"""


from dolfin import *
from utilities import *
import matplotlib.pyplot as plt

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

"""
three layer problem with a single production well

"""

# Create Mesh/Read Mesh from file
# mesh = UnitSquareMesh(10,10)
mesh_file_name = "grids/unitsqtriangle/unitsq_tri"
mesh = Mesh(mesh_file_name +".xml")
subdomains = MeshFunction("size_t", mesh, mesh_file_name+"_physical_region.xml")
bndry  = MeshFunction("size_t", mesh, mesh_file_name+"_facet_region.xml")
n = FacetNormal(mesh)

num_cell = mesh.num_cells()

KSpace = FunctionSpace(mesh,'DG', 0)
phi = Function(KSpace)
phi.vector()[:] = gaussianField((0.1,0.4), num_cell)
phi.rename("phi","poro")


k_ = Function(KSpace)
#k_.vector()[:] = 1e-8/(1000*9.81)
k_.vector()[:] = perm(phi.vector()[:])[:]
k_.rename("k", "perm")


# Save phi and k to the file
File("phi.pvd") << phi
File("k.pvd")  << k_


# # Define permeability expression and matrix
# phi = MeshFunction("double", mesh, "data/phi.xml.gz")
# p = phi.array()[:]
# kxx = MeshFunction("double", mesh, "data/kxx.xml.gz")
# kyy = MeshFunction("double", mesh, "data/kyy.xml.gz")
# tau_ = 0.81
# dp = 10e-6
# perm =  p**3*(dp)**2/(tau_*72*(1.- p)**2)
# kxx.array()[:] = perm[:]
# kyy.array()[:] = perm[:]


#from cpp_code import kcppcode
# from cpp_code import cppcode2
#k_ = CompiledExpression(compile_cpp_code(kcppcode).K3(),
#                       c00=kxx, c01=kyy, degree=0)

#phi = k_[0]
#Kinv = as_matrix(((1./k_[0], 0), (0, 1./k_[1])))
#Kinv = 1/k_[0]


# Parameters
day= 86400
MPa = 1e6
E = 100*MPa     #[Pa] Young's modulus
K  = 1e-8       # [m^2] Permeability
Kf = 2.0e9      # [Pa]
Ks = 1e10       # [Pa]
      
# phi = 0.375
alpha = 1.0     # Biot Coefficient
nu = 0.35 

M = Mb(phi,Kf,Ks) # [Pa] Biot Modulus

#mu = 1.002e-3       # [Pa s]
PL = -1.0e4      # [Pa] Vertical load
dt = dt = 0.005*day        # [s]

#K = Constant(K/(1000*9.81))
#K = Constant(K/mu)
#K = Constant(10e-5)


# s0 = Constant(1.0/M)
s0 = (1.0/M)
#s0 = MeshFunction("double", mesh, mesh.topology().dim(),0.0)
#s0.vector()[:] = (1.0/M)



# Class representing the intial conditions
class kxx(UserExpression):
    def __init__(self, kxx, **kwargs):
        #random.seed(2 + MPI.rank(MPI.comm_world))
        super().__init__(**kwargs)
        self.kxx = kxx
        

    def eval_cell(self, values, x, ufc_cell):
        index = ufc_cell.index
        values[0] = 1./self.kxx[index]  # s0
        
    def value_shape(self):
        return ()


# kxx = MeshFunction("double", mesh, mesh.topology().dim(),0.0)
        
    

# # Class representing the intial conditions
# class s01(UserExpression):
#     def __init__(self, phi, **kwargs):
#         #random.seed(2 + MPI.rank(MPI.comm_world))
#         super().__init__(**kwargs)
#         self.phi = phi
#         self.M = Mb(phi,Kf,Ks)
#         self.s0 = 1.0/M
        

#     def eval_cell(self, values, x, ufc_cell):
#         index = ufc_cell.index
#         values[0] = self.s0[index]  # s0
        
#     def value_shape(self):
#         return ()


#s0 = s0(p)
#Kinv = kxx(perm)
        
    

G = Constant(40*MPa)
lmbda = Constant(40*MPa)

#G = Constant(E/(2.0*(1.0+nu)))
#lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))

#lmbda = 40*MPa  # Lame's coefficient [Pa]
#G = 40*MPa  # \mu Lame's coefficient [Pa

# Class for interfacing with the Newton solver
class Biots3P(NonlinearProblem):
    def __init__(self, a, L, bc):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
        self.bc = bc

    def F(self, b, x):
        assemble(self.L, tensor=b)
        for bc in self.bc:
            bc.apply(b)

    def J(self, A, x):
        assemble(self.a, tensor=A)
        for bc in self.bc:
            bc.apply(A)



# Class representing the intial conditions
class InitialConditions(UserExpression):
    def __init__(self, **kwargs):
        #random.seed(2 + MPI.rank(MPI.comm_world))
        super().__init__(**kwargs)

    def eval(self, values, x):
        values[0] = 0.0  # ux
        values[1] = 0.0  # uy
        values[2] = 0.0  # qy
        values[3] = 0.0  # qy
        values[4] = 500.0  # p
        
    def value_shape(self):
        return (5,)
    


#kappa = K(materials, k_0, k_1, degree=1)
#print(kappa(0.5,0.5))

# K = K(subdomains=subdomains,degree=1)

#phi = phi(subdomains=subdomains,degree=1)

K = Constant(1e-8/(1000*9.81))
#Kinv = inv(k_)

Kinv = as_matrix(((1./k_, 0), (0, 1./k_)))

#phi = Constant(0.35)

""" Interpolate the K and Phi to check the interpolation 
U = FunctionSpace(mesh,"Lagrange",1)
#sol = Function(U)
#sol.interpolate(perm)
#filek = File("K.pvd")
#filek << sol
"""

I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I


# build function space
U = VectorElement("CG", mesh.ufl_cell(),2, 2) # displacement
V = FiniteElement("RT", mesh.ufl_cell(),1,2) # velocity
#V = FiniteElement("BDM", mesh.ufl_cell(),1) # velocity
P = FiniteElement("DG", mesh.ufl_cell(),0) # pressure
W = FunctionSpace(mesh, MixedElement([U,V,P]))

# Initial Condition
w0 = Function(W)  # previous time step solution
(u0, q0, p0) = split(w0)

# Define trail and test functions
dw = TrialFunction(W)
w = Function(W)  # current time step solution
(u,q,p) = split(w)
(v,z,y) = TestFunctions(W)

# Create intial conditions and interpolate
w_init = InitialConditions(degree=2)
w0.interpolate(w_init)
w.interpolate(w_init)


# Define bondary conditions

"""

left = 1
right = 2
top = 4
bottom = 3

"""

def getBCs():
    bcs = []
    
    # flux and displacement    
    # Left boundary
    bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,1)) # No normal displacement for solid on left
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0), bndry, 1)) 

    
    # Right boundary
    bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,2)) # No normal displacement for solid on right
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0), bndry,2))
    
    # top
    #ux_top = DirichletBC(W.sub(0).sub(0), 0.0, bndry,26) # No normal displacement for solid on left
    #uy_top = DirichletBC(W.sub(0).sub(1), 0.0, bndry,26) # No normal displacement for solid on bottom
    # bcs.append(DirichletBC(W.sub(1), (0.0,0.0), bndry,4)) 
    #p_top = DirichletBC(W.sub(2), 0.0, bndry,26)
    
    # bottom boundary
    bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, bndry,3)) # No normal displacement for solid on bottom
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0),  bndry,3))

    return bcs


# Volume force/source
f = Constant((0.0,0.0))
# Load/traction
t = Constant((0.0,PL))
g = Constant(0.0)

theta = 1.0
umid = (1.0-theta)*u0 + theta*u
pmid = (1.0-theta)*p0 + theta*p
qmid = (1.0-theta)*q0 + theta*q

k = Constant(dt)
dpdt = (p-p0)/k
dudt = (u-u0)/k

pbc = Constant(20.0)

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
dss = Measure('ds', domain=mesh, subdomain_data=bndry)
#dss = ds(subdomain_data = bndry)
# Weak form for elastic 
L0 = inner(sigma(umid), sym(grad(v)))*dx  -  alpha*pmid*div(v)*dx \
    - inner(f, v)*dx - inner(t, v)*dss(4) + alpha*inner(pbc*n,v)*dss(1)

# Darcy velocity
L1 =(inner(Kinv*qmid, z) - pmid*div(z))*dx - inner(z, pbc*n)*dss(1)

# Weak form for mass conservation
L2 = (s0*dpdt + alpha*div(dudt) + div(qmid) -g)*y*dx 
F = L0 + L1 + L2 

J = derivative(F,w, dw)


problem = Biots3P(J,F, getBCs())
solver = NewtonSolver()
solver.parameters["linear_solver"] = "gmres"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6
solver.parameters["maximum_iterations"] = 50


#solver.parameters["preconditioner"] = "ilu"


# Output file
ufile = XDMFFile("results/numerical/u.xdmf")
pfile = XDMFFile("results/numerical/p.xdmf")

# uFile = File("results/numerical/u.pvd")
# pFile = File("results/numerical/p.pvd")

u_, v_,p_ = w.split(deepcopy=True)
u_.rename("u","displacement")
p_.rename("p","pressure")
ufile.write(u_,0.0)
pfile.write(p_,0.0)

# uFile << (u_,0.0)
# pFile << (p_,0.0)


# Step in time
t = 0.0
T = 0.5*day

while (t < T):

    t += dt
    # Solve it
    solver.solve(problem, w.vector())

    # Save solution to file (VTK)
#    if t%8640==0:
    u_, v_,p_ = w.split(deepcopy=True)
    u_.rename("u","displacement")
    p_.rename("p","pressure")
    ufile.write(u_,t)
    pfile.write(p_,t)

    # Update previous solution
    w0.assign(w)
    

ufile.close()
pfile.close()

