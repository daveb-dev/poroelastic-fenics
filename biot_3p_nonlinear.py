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
        values[4] = 0.0

    def value_shape(self):
        return (5,)


# Create mesh
L = 10.0
nx = 10
ny = 60
mesh = RectangleMesh.create([Point(0.0,0.0), Point(1.0, L)],[nx,ny], CellType.Type.triangle)
n = FacetNormal(mesh)

# Parameters
MPa = 1e6
#E = 100*MPa        #[Pa] Young's modulus
K  = 1e-8           # [m^2] Permeability
Kf = 2.0e9      # [Pa]
cf = 1./Kf        
phi = 0.375
alpha = 1.0      # Biot Coefficient
nu = 0.35 



M = Mc(phi,cf) # [Pa] Biot Modulus

mu = 1.002e-3       # [Pa s]

L = 10.0            # [m]
PL = -1e5        # [Pa] Vertical load
dt = 1e-1        # [s]



K = Constant(K/(1000*9.81))
K = Constant(K/mu)
#K = Constant(10e-5)


s0 = Constant(1.0/M)

G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))


# Define Dirichlet boundary
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], L)

top = top()

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
top.mark(boundary_subdomains,1)


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
# displacement
bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "near(x[0], 0.0) || near(x[0], 1.0)") # No normal displacement for solid on sides
bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)")                    # No normal displacement for solid on bottom
#flux boundary
bcf = DirichletBC(W.sub(1), (0.0,0.0), "near(x[1], 0.0) || near(x[0], 0.0)|| near(x[0], 1.0)")   # Zero fluid pressure (drained) on top, natural bc
# Pressure
bcp = DirichletBC(W.sub(2), 0.0, "near(x[1], 10.0)")   # Zero fluid pressure (drained) on top, natural bc
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
solver.parameters["relative_tolerance"] = 1e-4


# Output file
filev = File("biot3p-nonlinear/u.pvd", "compressed")
filep = File("biot3p-nonlinear/p.pvd", "compressed")
filev<<(w.split()[0],0.0)
filep<<(w.split()[2],0.0)

# Step in time
t = 0.0
T = 12*dt;
count = 0
while (t < T):

    t += dt
    # Solve it
    #solve(a == L, w, bcs)
    solver.solve(problem, w.vector())

    # Save solution to file (VTK)
    if count%5==0:
        filev<<(w.split()[0],t)
        filep<<(w.split()[2],t)

    
    # Update previous solution
    w0.assign(w)
    count = count + 1

