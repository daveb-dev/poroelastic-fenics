from dolfin import *
from utilities import *
import matplotlib.pyplot as plt

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True


"""
three layer problem with a single production well

"""


# Parameters
MPa = 1e6
day = 60*60*24
E = 100*MPa     #[Pa] Young's modulus
K  = 1e-8       # [m^2] Permeability
Kf = 2.0e9      # [Pa]
Ks = 1e10       # [Pa]
      
phi = 0.375
alpha = 1.0     # Biot Coefficient
nu = 0.35 

M = Mb(phi,Kf,Ks) # [Pa] Biot Modulus
#mu = 1.002e-3       # [Pa s]
PL = -20*MPa      # [Pa] Vertical load

dt = 60  #0.01*day      # [s]

#K = Constant(K/(1000*9.81))
#K = Constant(K/mu)
#K = Constant(10e-5)

s0 = Constant(1.0/M)

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
        values[2] = 0.0  # qx
        values[3] = 0.0  # qy
        values[4] = 20*MPa  # p
        
    def value_shape(self):
        return (5,)




""" TODO: Implement K for layers """

class K(UserExpression):
    def __init__(self, subdomains, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.k0 = 1e-8/(1000*9.81)
        self.k1 = 8e-8/(1000*9.81)
        self.k2 = 1e-8/(1000*9.81)

    def eval_cell(self,values, x, cell):
        if self.subdomains[cell.index] == 6: # top subdomain
            values[0] = self.k0
        elif self.subdomains[cell.index] == 7: # reservoir subdomain
            values[0] = self.k1
        else:                              # bottom subdomain
            values[0] = self.k2
            
    def value_shape(self):
        return ()
            
            
class phi(UserExpression):
    def __init__(self, subdomains, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.phi0 = 0.1
        self.phi1 = 0.35
        self.phi2 = 0.1

    def eval_cell(self,values, x, cell):
        if self.subdomains[cell.index] == 6: # top subdomain
            values[0] = self.phi0
        elif self.subdomains[cell.index] == 7: # reservoir subdomain
            values[0] = self.phi1
        else:                              # bottom subdomain
            values[0] = self.phi2
                
    def value_shape(self):
        return ()


mesh_file_name = "grid/example01"
mesh = Mesh(mesh_file_name +".xml")
subdomains = MeshFunction("size_t", mesh, mesh_file_name+"_physical_region.xml")
bndry  = MeshFunction("size_t", mesh, mesh_file_name+"_facet_region.xml")
n = FacetNormal(mesh)

dx = dx(subdomain_data=subdomains)

#dx = Measure('dx')[subdomains]
#dss = Measure('ds')[bndry]
dss = ds(subdomain_data=bndry)


#kappa = K(materials, k_0, k_1, degree=1)
#print(kappa(0.5,0.5))

K = K(subdomains=subdomains,degree=1)
phi = phi(subdomains=subdomains,degree=1)



""" Interpolate the K and Phi to check the interpolation 
U = FunctionSpace(mesh,"Lagrange",1)
#sol = Function(U)
#sol.interpolate(perm)
#filek = File("K.pvd")
#filek << sol

"""


## Define Dirichlet boundary
#class top(SubDomain):
#    def inside(self, x, on_boundary):
#        return on_boundary and near(x[1], b)
#
#class bottom(SubDomain):
#    def inside(self, x, on_boundary):
#        return on_boundary and near(x[1], 0)
#
#top = top()
#bottom = bottom()
#
## Create mesh function over the cell facets
#boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
#boundary_subdomains.set_all(0)
#top.mark(boundary_subdomains,1)
##bottom.mark(boundary_subdomains,1)
#

I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I


# build function space
U = VectorElement("CG", mesh.ufl_cell(),2, 2) # displacement
V = FiniteElement("RT", mesh.ufl_cell(),1) # velocity
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

# Well
bcu1 = DirichletBC(W.sub(0).sub(0), 0.0, bndry,1) # No normal displacement for solid on left
#bcu1 = DirichletBC(W.sub(0).sub(1), 0.0, bndry,1) # No normal displacement for solid on bottom
#bcf1 = DirichletBC(W.sub(1), (0.0,0.0), bndry,1 ) # flux
#bcp1 = DirichletBC(W.sub(2), 0.0, bndry, 1 )  # pressure 


# Left
bcu5 = DirichletBC(W.sub(0).sub(0), 0.0, bndry,5) # No normal displacement for solid on left
#bcu5 = DirichletBC(W.sub(0).sub(1), 0.0, bndry,5) # No normal displacement for solid on bottom
bcf5 = DirichletBC(W.sub(1), (0.0,0.0), bndry,5) 
#bcp5 = DirichletBC(W.sub(2), 0.0, bndry, 5 )


# Right
bcu3 = DirichletBC(W.sub(0).sub(0), 0.0, bndry,3) # No normal displacement for solid on left
#bcu3 = DirichletBC(W.sub(0).sub(1), 0.0, bndry,3) # No normal displacement for solid on bottom
bcf3 = DirichletBC(W.sub(1), (0.0,0.0), bndry,3) 
#bcp3 = DirichletBC(W.sub(2), 0.0, bndry,3 )


# top
#bcu4 = DirichletBC(W.sub(0).sub(0), 0.0, bndry,4) # No normal displacement for solid on left
#bcu4 = DirichletBC(W.sub(0).sub(1), 0.0, bndry,4) # No normal displacement for solid on bottom
bcf4 = DirichletBC(W.sub(1), (0.0,0.0), bndry,4) 
#bcp4 = DirichletBC(W.sub(2), 0.0, bndry,4 )


# bottom
#bcu2 = DirichletBC(W.sub(0).sub(0), 0.0, bndry,2) # No normal displacement for solid on left
bcu2 = DirichletBC(W.sub(0).sub(1), 0.0,  bndry,2) # No normal displacement for solid on bottom
bcf2 = DirichletBC(W.sub(1), (0.0,0.0),  bndry,2) 
# bcp2 = DirichletBC(W.sub(2), 0.0, bndry,2 )

#bcs = [bcu1,bcp1, bcu5,bcf5,bcp5, bcu3,bcf3,bcp3, bcf4,bcp4, bcu2,bcf2,bcp2 ]
bcs = [bcu1, bcu5,bcf5, bcu3,bcf3, bcf4, bcu2,bcf2 ]
#
## displacement
#bc0 = DirichletBC(W.sub(0).sub(1), 0.0, bndry,2) # bottom No normal displacement for solid on left
#bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)") # No normal displacement for solid on bottom
#
##flux boundary
#bcf = DirichletBC(W.sub(1), (0.0,0.0), "near(x[1], 1.0) || near(x[1], 0.0) || near(x[0],0.0)")   # Zero fluid pressure (drained) on top & bottom, natural bc
##bcf = DirichletBC(W.sub(1), (0.0,0.0), "near(x[1], 1.0) || near(x[1], 0.0) || near(x[0],0.0)")   # Zero fluid pressure (drained) on top & bottom, natural bc
#
## Pressure
#bcp = DirichletBC(W.sub(2), 0.0, bndry, 1 )   # Zero fluid pressure (drained) on left and right, natural bc
##bcp2 = DirichletBC(W.sub(2), 0.0, "near(x[0], 1.0)" )   # Zero fluid pressure (drained) on left and right, natural bc
##bcs = [bc0, bcf, bcp]
#bcs = [bc0, bc1, bcf, bcp]

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

pbc = Constant(2.0*MPa) # BHP 2.0 MPa

#dss = ds(subdomain_data =bndry)
# Weak form for elastic 
L0 = inner(sigma(umid), sym(grad(v)))*dx  -  alpha*pmid*div(v)*dx \
    - inner(f, v)*dx - inner(t, v)*dss(4)  - alpha*inner(pbc*n,v)*dss(1)

# Darcy velocity
L1 =(inner(inv(K)*qmid, z) - pmid*div(z))*dx + inner(z, pbc*n)*dss(1)

# Weak form for mass conservation
L2 = (s0*dpdt + alpha*div(dudt) + div(qmid) -g)*y*dx 
F = L0 + L1 + L2 

J = derivative(F,w, dw)

problem = Biots3P(J,F, bcs)

#from IPython import embed; embed()

solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6


# Output file
ufile = XDMFFile("results/numerical/u.xdmf")
pfile = XDMFFile("results/numerical/p.xdmf")
u_, v_,p_ = w.split(deepcopy=True)
u_.rename("u","displacement")
p_.rename("p","pressure")
ufile.write(u_,0.0)
pfile.write(p_,0.0)


# Step in time
t = 0.0
T =  7200 #0.5*day

count = 1

while (t < T):

    t += dt
    # Solve it
    #solve(a == L, w, bcs)
    solver.solve(problem, w.vector())

    # Save solution to file (VTK)
    #if count%100==0:
    u_, v_,p_ = w.split(deepcopy=True)
    u_.rename("u","displacement")
    p_.rename("p","pressure")
    ufile.write(u_,t)
    pfile.write(p_,t)

    # Update previous solution
    w0.assign(w)
    count = count + 1
    
ufile.close()
pfile.close()

