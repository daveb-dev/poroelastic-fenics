from dolfin import *
from utilities import *
import matplotlib.pyplot as plt
#from IPython import embed; embed()

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

"""
three layer problem with a single production well

"""

# Parameters
MPa = 1e6
E = 100*MPa     #[Pa] Young's modulus
K  = 1e-8       # [m^2] Permeability
Kf = 2.0e9      # [Pa]
Ks = 1e10       # [Pa]
      
phi = 0.375
alpha = 1.0     # Biot Coefficient
nu = 0.35 

M = Mb(phi,Kf,Ks) # [Pa] Biot Modulus
#mu = 1.002e-3       # [Pa s]
PL = -1.0e4      # [Pa] Vertical load
dt =0.0001      # [s]

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
        values[2] = 0.0  # ux
        values[3] = 0.0  # qy
        values[4] = 0.0  # qy
        values[5] = 0.0  # qy
        values[6] = 1.0  # p
        
    def value_shape(self):
        return (7,)




""" TODO: Implement K for layers """

#class K(UserExpression):
#    def __init__(self, subdomains, **kwargs):
#        super().__init__(**kwargs)
#        self.subdomains = subdomains
#        self.k0 = 1e-8/(1000*9.81)
#        self.k1 = 3e-8/(1000*9.81)
#        self.k2 = 8e-8/(1000*9.81)
#
#    def eval_cell(self,values, x, cell):
#        if self.subdomains[cell.index] == 6: # top subdomain
#            values[0] = self.k0
#        elif self.subdomains[cell.index] == 7: # reservoir subdomain
#            values[0] = self.k1
#        else:                              # bottom subdomain
#            values[0] = self.k2
#            
#            
#class phi(UserExpression):
#    def __init__(self, subdomains, **kwargs):
#        super().__init__(**kwargs)
#        self.subdomains = subdomains
#        self.phi0 = 0.1
#        self.phi1 = 0.35
#        self.phi2 = 0.1
#
#    def eval_cell(self,values, x, cell):
#        if self.subdomains[cell.index] == 6: # top subdomain
#            values[0] = self.phi0
#        elif self.subdomains[cell.index] == 7: # reservoir subdomain
#            values[0] = self.phi1
#        else:                              # bottom subdomain
#            values[0] = self.phi2
#
#



mesh_file_name = "grid/example02a"
mesh = Mesh(mesh_file_name +".xml")
subdomains = MeshFunction("size_t", mesh, mesh_file_name+"_physical_region.xml")
bndry  = MeshFunction("size_t", mesh, mesh_file_name+"_facet_region.xml")


n = FacetNormal(mesh)


#kappa = K(materials, k_0, k_1, degree=1)
#print(kappa(0.5,0.5))

#K = K(subdomains=subdomains,degree=1)
#phi = phi(subdomains=subdomains,degree=1)

K = Constant(1e-8/(1000*9.81))
phi = Constant(0.35)

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
U = VectorElement("CG", mesh.ufl_cell(),2, 3) # displacement
V = FiniteElement("RT", mesh.ufl_cell(),1,3) # velocity
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
well = 0
top = 26
bottom = 29
front = 30
back = 27
left = 28
right = 31

"""

def getBCs():
    bcs = []
    # flux and displacement    
    # Left boundary
    bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,28)) # No normal displacement for solid on left
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0), bndry,28)) 

    
    # Right boundary
    bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,31)) # No normal displacement for solid on right
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0), bndry,31))
    
    # front boundary
    bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, bndry,30)) # No normal displacement for solid on front
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0), bndry,30)) 

    
    # back boundary
    bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, bndry,27)) # No normal displacement for solid on back
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0), bndry,27))

    
    # top
    #ux_top = DirichletBC(W.sub(0).sub(0), 0.0, bndry,26) # No normal displacement for solid on left
    #uy_top = DirichletBC(W.sub(0).sub(1), 0.0, bndry,26) # No normal displacement for solid on bottom
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0), bndry,26)) 
    #p_top = DirichletBC(W.sub(2), 0.0, bndry,26)
    
    # bottom boundary
    bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,29)) # No normal displacement for solid on bottom
    bcs.append(DirichletBC(W.sub(1), (0.0,0.0,0.0),  bndry,29))

    # Well
    #bcs.append(DirichletBC(W.sub(0).sub(0), 0.0, bndry,0)) # No normal displacement for solid on left
    #bcs.append(DirichletBC(W.sub(0).sub(1), 0.0, bndry,0)) # No normal displacement for solid on left

    return bcs


# Volume force/source
f = Constant((0.0,0.0,0.0))
# Load/traction
t = Constant((0.0,0.0,PL))
g = Constant(0.0)

theta = 1.0
umid = (1.0-theta)*u0 + theta*u
pmid = (1.0-theta)*p0 + theta*p
qmid = (1.0-theta)*q0 + theta*q

k = Constant(dt)
dpdt = (p-p0)/k
dudt = (u-u0)/k

pbc = Constant(5000.0)

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
dss = Measure('ds', domain=mesh, subdomain_data=bndry)
#dss = ds(subdomain_data = bndry)
# Weak form for elastic 
L0 = inner(sigma(umid), sym(grad(v)))*dx  -  alpha*pmid*div(v)*dx \
    - inner(f, v)*dx - inner(t, v)*dss(26) + alpha*inner(pbc*n,v)*dss(0)

# Darcy velocity
L1 =(inner(inv(K)*qmid, z) - pmid*div(z))*dx - inner(z, pbc*n)*dss(0)

# Weak form for mass conservation
L2 = (s0*dpdt + alpha*div(dudt) + div(qmid) -g)*y*dx 
F = L0 + L1 + L2 

J = derivative(F,w, dw)


problem = Biots3P(J,F, getBCs())
solver = NewtonSolver()
solver.parameters["linear_solver"] = "gmres"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-3
solver.parameters["maximum_iterations"] = 50
solver.parameters["preconditioner"] = "ilu"


# Output file
ufile = XDMFFile("results/numerical/u.xdmf")
pfile = XDMFFile("results/numerical/p.xdmf")

uFile = File("results/u.pvd")
pFile = File("results/p.pvd")

u_, v_,p_ = w.split(deepcopy=True)
u_.rename("u","displacement")
p_.rename("p","pressure")
ufile.write(u_,0.0)
pfile.write(p_,0.0)

uFile << (u_,0.0)
pFile << (p_,0.0)


# Step in time
t = 0.0
T = 2
count = 1
while (t < T):

    t += dt
    # Solve it
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

