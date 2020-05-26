from dolfin import *
import matplotlib.pyplot as plt

# Create mesh
nx = 10
ny = 60
mesh = RectangleMesh.create([Point(0.0,0.0), Point(2.0,15.0)],[nx,ny], CellType.Type.quadrilateral)

# Parameters
E = 100e6         #[Pa] Young's modulus
alpha = 0.79      # Biot Coefficient
nu = 0.25
M     = 12.0e9   # [Pa] Biot Modulus
k  = 1.875e-5        # [m^2] Permeability
mu = 1 #1.002e-3       # [Pa s]
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
print("s0", float(s0))

mu = 1.0e-5 # [m^2/s]
G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))

# Mass density
rho = Constant(2600)


# Define permeability k aka hydraulic conductance
#epsilon = 1.e-8
#k = Expression("(x[1] > 0.25 && x[1] <= 0.75) ? epsilon : 1.0",
#               epsilon=epsilon, degree=1)

#k = Constant(1.0e-14) # m^2

#plot(k, interactive=True, mesh=mesh)

# Define other material parameters
#s0 = Constant(s0)
#alpha = Constant(alpha)
#dt = Constant(1.0)


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




n = FacetNormal(mesh)
#Fx = Expression("near(x[1], 15.0) ? -1.0 : 0.0", degree=1)
#Fx = Constant(-1000)

# build function space
V = VectorElement("CG", mesh.ufl_cell(),2) # displacement
Q = FiniteElement("CG", mesh.ufl_cell(),1) # pressure

M = MixedElement([V,Q])
W = FunctionSpace(mesh, M)


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
w.interpolate(w_init)

# Volume force/source
f = Constant((0.0,-PL))
s = Constant(0.0)

rho_g = Constant((0.0, -rho*9.81))
rho_g = Constant((0.0, 0.0))

dss = ds(subdomain_data =boundary_subdomains)
# Weak form for elastic 
L0 = inner(sigma(u) - alpha*p*I, sym(grad(v)))*dx \
    - inner(f, v)*dss(1) - inner(rho_g, v)*dx \
    #- inner(rho_g, v)*dx \

# Weak form for mass conservation
L1 = (s0*p + alpha*div(u))*q*dx \
    + dt*inner(k*grad(p), grad(q))*dx \
    - (s0*p0 + alpha*div(u0) + dt*s)*q*dx

F = L0 + L1 

(a, L) = system(F)

# Define bondary conditions
# displacement
bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "near(x[0], 0.0) || near(x[0], 2.0)") # No normal displacement for solid on sides
bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)")                    # No normal displacement for solid on bottom
# Pressure
bcp = DirichletBC(W.sub(1), 0.0, "near(x[1], 15.0)")                           # Zero fluid pressure (drained) on top, natural bc
bcs = [bc0, bc1, bcp]


# Output file
filev = File("2phase/u.pvd", "compressed")
filep = File("2phase/p.pvd", "compressed")

# Step in time
t = 0.0
T = 1000*dt;
count = 0
while (t < T):

    if count%50==0:
        filev<<(w.split()[0],t)
        filep<<(w.split()[1],t)

    # Solve it
    solve(a == L, w, bcs)
    #solve(F, w, bcs)

    
    w0.assign(w)
    t += dt
    count = count + 1

#    # Plot fields
#    (u, p) = w.split(deepcopy=True)
#    fig1 = plt.figure()
#    plot(u, title="u")
#    fig2 = plt.figure()
#    cs = plot(p, title="p")
#    cbar = fig2.colorbar(cs)
#    plt.show()
#
##plot_over_line(p, k)
#
#
#plt.rcParams['figure.figsize'] = 20, 15
#plt.rcParams['lines.linewidth'] = 2
#plt.rcParams['image.cmap'] = 'jet'

def plot_over_line(p,k):
    import numpy as np
    plt.rcParams['figure.figsize'] = 20, 15
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['image.cmap'] = 'jet'
    points = np.linspace(0,1,100)
    values = [p(0.75,z) for z in points]
    plt.plot(points, values, 'b-', label="p")
    #values = [k(0.75,z) for z in points]
    #plt.plot(points, values, 'r--', label="k")
    plt.grid(True)
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("p")
    plt.show()


