from dolfin import *
import matplotlib.pyplot as plt

# Create mesh
nx = 10
ny = 60
mesh = RectangleMesh.create([Point(0.0,0.0), Point(1.0,15.0)],[nx,ny], CellType.Type.quadrilateral)
n = FacetNormal(mesh)

# Parameters
E = 100e6         #[Pa] Young's modulus
alpha = 1.0      # Biot Coefficient
nu = 0.25
M     = 12.0e9   # [Pa] Biot Modulus
k  = 1.e-10        # [m^2] Permeability
mu = 1.002e-3       # [Pa s]

k = k/mu
L = 15.0            # [m]
PL = -1e3        # [Pa] Vertical load
dt = 1e-1        # [s]
rho = 1000 
G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))
phi = 0.375
cf = -1.02
cs = 0.02
c0 = phi*cf + (alpha - phi)*cs
s0 = Constant(1.0/M)
s0 = Constant(4.5e-10)
print("s0", float(s0))

mu = 1.0e-5 # [m^2/s]
G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))

# Mass density
rho = Constant(2600)

# Initial Condition
class InitialConditions(UserExpression):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)

    def eval(self, values, x):
        values[0] = 0           # x displacement
        values[1] = 0           # y displacement
        values[2] = 0           #pressure

    def value_shape(self):
        return (3,)


# Define Dirichlet boundary
class right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 2)


class left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)


class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 15.0)


class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)

right = right()
left = left()
top = top()
bottom = bottom()

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
top.mark(boundary_subdomains,1)


I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I

# build function space
V = VectorElement("CG", mesh.ufl_cell(),2) # displacement
Q = FiniteElement("CG", mesh.ufl_cell(),1) # pressure

#M = MixedElement([V,Q])
W = FunctionSpace(mesh, V*Q)

# Initial Condition
w0 = Function(W)
(u0,p0) = split(w0)
w_init = InitialConditions(degree=1)
w0.interpolate(w_init)


# Define trail and test functions
(u,p) = TrialFunctions(W)
(v,q) = TestFunctions(W)

# Variational formulation
w = Function(W)
#w.interpolate(w_init)

# Volume force/source
f = Constant((0.0,0.0))
t = Constant((0.0,-PL))
g = Constant(0.0)

dss = ds(subdomain_data =boundary_subdomains)

k = Constant(dt)
dpdt = (p - p0)/k
dudt = (u - u0)/k
theta = Constant(0.5)
um = theta*u + (1 - theta)*u0
pm = theta*p + (1 - theta)*p0


# Weak form for elastic 
#L0 = (inner(sigma(um),grad(v)) + alpha*pm*div(v) - inner(f, v))*dx - inner(t, v)*dss(1) 

# Weak form for mass conservation
#L1 = (s0*dpdt*q  + alpha*div(dudt)*q + inner(k*grad(pm), grad(q)) - g*q)*dx   

# term k/mu grad(p)*q*ds will be zero here

#L0 = (inner(sigma(um), sym(grad(v))) + alpha * dot(grad(pm), v) - dot(f,v))*dx - inner(t,v)*dss(1)
#L1 = (s0*dpdt*q  + alpha*div(dudt)*q + inner(k*grad(pm), grad(q)) - g*q)*dx   


#F = L0 + L1 

# Variational formulation
F = (inner(sigma(um), sym(grad(v))) - alpha*pm*div(v) - inner(f, v)
     - s0*dpdt*q - alpha*div(dudt)*q - inner(k*grad(pm), grad(q)) - g*q)*dx \
     - (inner(t, v))*dss(1)

(a, L) = system(F)

# Define bondary conditions
# displacement
bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "near(x[0], 0.0) || near(x[0], 1.0)") # No normal displacement for solid on sides
bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)")                    # No normal displacement for solid on bottom
# Pressure
bcp = DirichletBC(W.sub(1), 0.0, "near(x[1], 15.0)")                           # Zero fluid pressure (drained) on top, natural bc
bcs = [bc0, bc1, bcp]


# Output file
filev = File("2phase/u.pvd", "compressed")
filep = File("2phase/p.pvd", "compressed")

# Step in time
t = 0.0
T = 10*dt;
count = 0
while (t < T):

    if count%1==0:
        filev<<(w.split()[0],t)
        filep<<(w.split()[1],t)

    # Solve it
    solve(a == L, w, bcs)

    w0.assign(w)
    t += dt
    count = count + 1
