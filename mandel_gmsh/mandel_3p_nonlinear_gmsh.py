from dolfin import *
from utilities import *
import matplotlib.pyplot as plt

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

    
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
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0
        values[3] = 0.0
        values[4] = 0.0 #4950.495049504951
    def value_shape(self):
        return (5,)


# Create mesh
a = 1.0
b = 1.0
nx = 30
ny = 30
#mesh = RectangleMesh.create([Point(0.0,0.0), Point(a, b)],[nx,ny], CellType.Type.triangle, "right/left")
path = "grids/case3/"
mesh = Mesh(path+"mesh.xml")
n = FacetNormal(mesh)

# Parameters
MPa = 1e6
E = 100*MPa     #[Pa] Young's modulus
K  = 1e-8       # [m^2] Permeability
Kf = 2.0e9      # [Pa]
cf = 1./Kf        
phi = 0.3
alpha = 1.0     # Biot Coefficient
nu = 0.25 



M = Mc(phi,cf) # [Pa] Biot Modulus
mu = 1.002e-3       # [Pa s]
PL = -1.0e4      # [Pa] Vertical load
dt = 1       # [s]



K = Constant(K/(1000*9.81))
#K = Constant(K/mu)
#K = Constant(10e-5)


s0 = Constant(1.0/M)

G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))



#lmbda = 40*MPa  # Lame's coefficient [Pa]
#G = 40*MPa  # \mu Lame's coefficient [Pa


# Define Dirichlet boundary
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], b)

class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0)

top = top()
bottom = bottom()

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
#boundary_subdomains = MeshFunction("size_t", mesh, path+'mesh_physical_region.xml')
#boundary_subdomains.set_all(0)
#top.mark(boundary_subdomains,1)
#bottom.mark(boundary_subdomains,1)




I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I


# build function space
cell = mesh.ufl_cell()
U = VectorElement("CG", cell ,2, 2) # displacement
V = FiniteElement("RT", cell ,1,2) # velocity
#V = FiniteElement("BDM", mesh.ufl_cell(),1) # velocity
P = FiniteElement("DG", cell ,0,2) # pressure
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
# displacement
#bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "near(x[0], 0.0)") # No normal displacement for solid on left
bc0 = DirichletBC(W.sub(0).sub(0), 0.0, boundary_subdomains, 1) 
#bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)") # No normal displacement for solid on bottom
bc1 = DirichletBC(W.sub(0).sub(1), 0.0, boundary_subdomains, 2) # No normal displacement for solid on bottom

#flux boundary
#bcf = DirichletBC(W.sub(1), (0.0,0.0), "near(x[1], 1.0) || near(x[1], 0.0) || near(x[0],0.0)")   # Zero fluid pressure (drained) on top & bottom, natural bc
bcf = DirichletBC(W.sub(1), (0.0,0.0), boundary_subdomains,3)   # no flux 

# Pressure
#bcp = DirichletBC(W.sub(2), 0.0, "near(x[0], 1.0)" )   # Zero fluid pressure (drained) on left and right, natural bc
bcp = DirichletBC(W.sub(2), 0.0, boundary_subdomains,3 )   # Zero fluid pressure (drained) on left and right, natural bc
#bcs = [bc0, bcf, bcp]
bcs = [bc0, bc1, bcf, bcp]

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

dss = ds(subdomain_data =boundary_subdomains)
# Weak form for elastic 
L0 = inner(sigma(umid), sym(grad(v)))*dx  -  alpha*pmid*div(v)*dx \
    - inner(f, v)*dx - inner(t, v)*dss(1) 

# Darcy velocity
L1 =(inner(inv(K)*qmid, z) - pmid*div(z))*dx

# Weak form for mass conservation
L2 = (s0*dpdt + alpha*div(dudt) + div(qmid) -g)*y*dx 
F = L0 + L1 + L2 

J = derivative(F,w, dw)

problem = Biots3P(J,F, bcs)
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6


# Output file
filev = File("mandel3p-nonlinear/u.pvd", "compressed")
filep = File("mandel3p-nonlinear/p.pvd", "compressed")
u_, v_,p_ = w.split(deepcopy=True)
p0_ =  p_.vector().max()
print("P0 = ", p0_)
#p_.vector()[:] = p_.vector()[:]/p0_
filev<<(w.split()[0],0.0)
p_.vector()[:] = p_.vector()[:]/4950.495049504951
filep<<(p_,0.0)

# Step in time
t = 0.0
T = 705
count = 1
while (t < T):

    t += dt
    # Solve it
    #solve(a == L, w, bcs)
    solver.solve(problem, w.vector())

    # Save solution to file (VTK)
    if count%10==0:
        u_, v_,p_ = w.split(deepcopy=True)
        p_.vector()[:] = p_.vector()[:]/4950.495049504951
        filev<<(u_,t)
        filep << (p_, t)

    
    # Update previous solution
    w0.assign(w)
    count = count + 1

