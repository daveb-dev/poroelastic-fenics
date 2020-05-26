from dolfin import *
import matplotlib.pyplot as plt

# Create mesh
nx = 10
ny = 60
mesh = RectangleMesh.create([Point(0.0,0.0), Point(2.0,15.0)],[nx,ny], CellType.Type.triangle)
n = FacetNormal(mesh)

# Parameters
E = 100e6         #[Pa] Young's modulus
alpha = 0.79      # Biot Coefficient
nu = 0.25
M     = 12.0e9   # [Pa] Biot Modulus
K  = 1.875e-5        # [m^2] Permeability
mu = 1.002e-3       # [Pa s]
L = 15.0            # [m]
PL = -1e3        # [Pa] Vertical load
dt = 1e-1        # [s]


G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))
phi = 0.375
s0 = Constant(1.0/M)

mu = 1.0e-5 # [m^2/s]
G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))


# Define Dirichlet boundary
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 15.0)

top = top()

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
top.mark(boundary_subdomains,1)


I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I


# build function space
U = VectorElement("CG", mesh.ufl_cell(),1) # displacement
V = FiniteElement("RT", mesh.ufl_cell(),1) # velocity
#V = FiniteElement("BDM", mesh.ufl_cell(),1) # velocity
P = FiniteElement("DG", mesh.ufl_cell(),0) # pressure
W = FunctionSpace(mesh, MixedElement([U,V,P]))

# Initial Condition
w0 = Function(W)
(u0, q0, p0) = split(w0)

# Define trail and test functions
(u,q,p) = TrialFunctions(W)
(v,z,y) = TestFunctions(W)

# Variational formulation
w = Function(W)

# Volume force/source
f = Constant((0.0,0.0))
# Load/traction
t = Constant((0.0,-PL))
g = Constant(0.0)

theta = 0.5
umid = (1.0-theta)*u0 + theta*u
pmid = (1.0-theta)*p0 + theta*p
qmid = (1.0-theta)*q0 + theta*q

dpdt = (p-p0)/dt
dudt = (u-u0)/dt

dss = ds(subdomain_data =boundary_subdomains)
# Weak form for elastic 
L0 = inner(sigma(umid), sym(grad(v)))*dx  -  alpha*pmid*div(v)*dx \
    - inner(t, v)*dss(1) - inner(f, v)*dx 

# Darcy velocity
L1 =(inner(inv(K)*qmid, z) - pmid*div(z))*dx

# Weak form for mass conservation
L2 = (s0*dpdt + alpha*div(dudt) + div(qmid) -g)*y*dx 




F = L0 + L1 + L2 

(a, L) = system(F)

# Define bondary conditions
# displacement
bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "near(x[0], 0.0) || near(x[0], 2.0)") # No normal displacement for solid on sides
bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)")                    # No normal displacement for solid on bottom
#flux boundary
bcflux = DirichletBC(W.sub(1), Constant((0,0)), "on_boundary") 
# Pressure
bcp = DirichletBC(W.sub(2), 0.0, "near(x[1], 15.0)")                           # Zero fluid pressure (drained) on top, natural bc
bcs = [bc0, bc1, bcp]


# Output file
filev = File("3phase/u.pvd", "compressed")
filep = File("3phase/p.pvd", "compressed")

# Step in time
t = 0.0
T = 1000*dt;
count = 0
while (t < T):

    if count%50==0:
        filev<<(w.split()[0],t)
        filep<<(w.split()[2],t)

    # Solve it
    solve(a == L, w, bcs)
    #solve(F, w, bcs)

    
    w0.assign(w)
    t += dt
    count = count + 1



