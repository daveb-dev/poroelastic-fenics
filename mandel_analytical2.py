import numpy as np
from utilities import *
from fenics import *

# Parameters
MPa = 1e6
E = 100*MPa        #[Pa] Young's modulus
K  = 1e-8           # [m/s] Permeability
Kf = 2.0e9      # [Pa]
cf = 1./Kf        
phi = 0.3
alpha = 1.0      # Biot Coefficient
nu = 0.25 

lmbda = 40*MPa  # Lame's coefficient [Pa]
G = 40*MPa  # \mu Lame's coefficient [Pa



M = Mc(phi,cf) # [Pa] Biot Modulus
M = Mb(phi,Kf)
print("Mb=", M)
print("Mc=",Mc(phi,cf))
Ku = Ku(lmbda,G,M)
print("Ku=", Ku)
cm = CM(lmbda,G)
print("CM", cm)
c = c(K,M,cm) 
print("c",c)
#K = K/(1000*9.81)
mu = 1.002e-3       # [Pa s]
F = 1.0e4      # [Pa] Vertical load
dt = 1.0        # [s]


B = getB(M,Ku)
print("B=", B)


nu = 0.25
nu_u = nu_u(nu,B)

print("nu_u=", nu_u)


class AnalyticalSol(UserExpression):
    def __init__(self,F,c,B,G,a=1.0, b=1.0, nu=0.25,nu_u=0.5, alpha=1.0, t=0.0, **kwargs):
        super().__init__(**kwargs)
        self.F = F
        self.c = c
        self.t = t
        self.G = G
        self.nu = nu
        self.B = B
        self.nu_u = nu_u
        self.a = a
        self.b = b
        self.p0 = F*B *(1 + nu_u)/(3.0*a)
        self.ux_0_term = F*nu_u / (2.0*G*a)
        self.uy_0_term = -F*(1.-nu_u)/(2.0*G*a)
        self.tmp1 = F*nu/(2*G*a)
        self.tmp2 = F*nu_u/(G*a)
        self.tmp3 = F*(1-nu)/(2*G*a)
        self.tmp4 = F*(1-nu_u)/(G*a)
        #self.beta_i = [1.3276, 4.641544, 7.81184, 10.9655455, 14.113833, 17.259678]
        self.beta_i =getBeta()
        """ TODO: update beta_i """

    def eval(self,value, x):
        """ Evaluate the pressure expression """
        ux = 0.0
        ux_1 = 0.0
        ux_2 = 0.0
        uy = 0.0
        p = 0.0
        for beta_i in self.beta_i:
            beta_i = float(beta_i*np.pi/180)
            expo_term = exp(-((beta_i**2)*self.c*self.t/self.a**2)) 
            term1 = sin(beta_i)/(beta_i - sin(beta_i)*cos(beta_i))
            p =  term1 * (cos(beta_i*x[0]/self.a)-cos(beta_i)) * expo_term
            ux_1 =  term1* cos(beta_i) * expo_term
            ux_2 = (term1*cos(beta_i)/sin(beta_i))*sin(beta_i*x[0]/self.a)*expo_term

            uy = term1*cos(beta_i)* expo_term

        ux = (self.tmp1 - self.tmp2*ux_1)*x[0] + (F/G)*ux_2
        
        value[0] = (self.tmp1 - self.tmp2*ux_1)*x[0] + (F/G)*ux_2  # ux
        value[1] = (-self.tmp3 + self.tmp4*uy)*x[1]    # uy 
        value[2] = 2*self.p0*p        # pressure 
    
    def value_shape(self):
        return (3,)


a = 1.0
b = 1.0
nx=10
ny=10
mesh = RectangleMesh.create([Point(0.0,0.0), Point(a,b)],[nx,ny], CellType.Type.triangle) 

U = VectorElement("Lagrange", mesh.ufl_cell(),2, 2) # displacement 
P = FiniteElement("Lagrange", mesh.ufl_cell(),1) # pressure
W = FunctionSpace(mesh, U*P)

w = Function(W)

fileu = File("analytical_mandel/u_ana.pvd")
filep = File("analytical_mandel/p_ana.pvd")

sol = AnalyticalSol(F, c,B,G,a,b)
p0 = sol.p0
w.interpolate(sol)

u_, p_ = w.split()
p_.vector()[:] = p_.vector()[:]/p0
fileu << (u_,0)
filep << (p_,0)

t = 0
T = 1000*60
dt = 10*60
time_steps = [1,30,120,240,480]
for t_ in time_steps:
    # update
    t = t_*60
    sol.t = t
    
    w.interpolate(sol)
    u_, p_ = w.split(deepcopy=True)
    p_.vector()[:] = p_.vector()[:]/p0
    fileu << (u_,t)
    filep << (p_,t)

