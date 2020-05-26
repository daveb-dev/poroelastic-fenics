from dolfin import *
import matplotlib.pyplot as plt

# Create mesh
L = 70
nx = 10
mesh = IntervalMesh(nx,0.0,L)

# Parameters
E = 584.0        #[Pa] Young's modulus
alpha = Constant(1)      # Biot Coefficient
nu = 0.35
K  = 1.875e-5        # [m^2] Permeability
c = Constant(0.0)#10.0**(-10))
G = Constant(E/(2.0*(1.0+nu)))
lmbda = Constant(E*nu/((1.0+nu)*(1.0-2.0*nu)))


I = Identity(mesh.topology().dim())
def sigma(u):
    return 2.0*G*sym(grad(u)) + lmbda*div(u)*I




##Fx = Expression("near(x[1], 15.0) ? -1.0 : 0.0", degree=1)
##Fx = Constant(-1000)
#
# build function space
V = VectorElement("CG", mesh.ufl_cell(),2) # displacement
Q = FiniteElement("CG", mesh.ufl_cell(),1) # pressure

#M = MixedElement([V,Q])
W = FunctionSpace(mesh, V*Q)




# Define trail and test functions
(u,p) = TrialFunctions(W)
(v,q) = TestFunctions(W)

# Solutions at previous time step
w_ = Function(W)
u_, p_ = split(w_)
w = Function(W)  # current time step solution
# Initial Condition
#w0 = Function(W)
#(u0,p0) = split(w0)
#w_init = InitialConditions(degree=1)
#w0.interpolate(w_init)


k = Constant(0.5)
T = 10.0
dpdt = (p-p_)/k
dudt = (u-u_)/k

theta = 1.0  # theta=1 -> backward, theta=0.5 -> Crank-Nicolson 
um = theta*u + (1-theta)*u_
pm = theta*p + (1-theta)*p_
dim = mesh.topology().dim()

# Volume force/source
f = Constant((0.0,)*dim)
g = Constant(0.0)

# Essential boundary conditions
u_bar = Constant((0.0,)*dim)
p_bar = Expression("A*(1.0 - x[0]/L) + B", degree=1, A=13.0, B=0.0, L=L)

#p_bar = Constant(0.0)
n = FacetNormal(mesh)
t = -alpha*15*n
s = Constant(0.0)

# Variational formulation
F = (inner(sigma(um), grad(v)) + alpha*pm*div(v) - inner(f, v)
     + alpha*div(dudt)*q \
     - inner(K*grad(pm), grad(q)) \
     - c*dpdt*q \
     - g*q)*dx \
     - (inner(t, v) - s*q)*ds

(a, L) = system(F)



#rho_g = Constant((0.0, -rho*9.81))
#rho_g = Constant((0.0, 0.0))
#
#dss = ds(subdomain_data =boundary_subdomains)
## Weak form for elastic 
L0 = inner(sigma(u) - alpha*p*I, sym(grad(v)))*dx \
    - inner(t, v)*ds - inner(f, v)*dx 

# Weak form for mass conservation
L1 = (c*p + alpha*div(u))*q*dx \
    + dt*inner(k*grad(p), grad(q))*dx \
    - (c*p0 + alpha*div(u0) + dt*s)*q*dx

F = L0 + L1 

(a, L) = system(F)

# Define bondary conditions
#u_bar = Constant((0.0,)*dim)
# displacement
bc0 = DirichletBC(W.sub(0), u_bar, "near(x[0], 0.0) && on_boundary") # No normal displacement for solid on sides
#bc0 = DirichletBC(W.sub(0).sub(0), 0.0, "(near(x[0], 0.0) || near(x[0], 2.0)) && on_boundary") # No normal displacement for solid on sides
#bc1 = DirichletBC(W.sub(0).sub(1), 0.0, "near(x[1], 0.0)")                    # No normal displacement for solid on bottom
## Pressure
bcp = DirichletBC(W.sub(1), p_bar, "near(x[0], 70.0) && on_boundary")                           # Zero fluid pressure (drained) on top, natural bc
bcs = [bc0,bcp]


## Output file
fileu = File("1D2phase/u.pvd", "compressed")
filep = File("1D2phase/p.pvd", "compressed")
#filev = File("1D2phase/q.pvd", "compressed")

w = Function(W)

## Step in time
t = Constant(0.0) 
T = 1000*float(k);
count = 0
while (float(t)< T):

    if count%50==0:
        #filev<<(w.split()[0],t)
        #filep<<(w.split()[1],t)
        filep << w.split()[1]

    # Solve it
    solve(a == L, w, bcs)

    (u,p) = w.split(deepcopy=True)

    # update 
    t.assign(float(t) + k)
    w_.assign(w)
    count = count + 1

