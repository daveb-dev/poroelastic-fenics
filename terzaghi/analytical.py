import numpy as np
from utilities import *
from fenics import *


# cf = 5
L = 10.0
nx=10
ny=40
MPa = 1e6
E = 100*MPa
K = 1e-8  # [m/s]
Kf = 2.0e9  # [Pa]
Ks = 1e10   # [Pa] 
cf = 1./Kf
phi = 0.375
alpha = 1.0  # Biot coefficient
PL = 1.0e4  # use load in [Pa] 
lmbda = 40*MPa  # Lame's coefficient [Pa]
G = 40*MPa  # \mu Lame's coefficient [Pa
#nu = 0.35


print("Mc=",Mc(phi,cf))
M = Mb(phi,Kf)
#M = Mc(phi,cf)
print("Mb=", M)
Ku = Ku(lmbda,G,M)
print("Ku=", Ku)
cm = CM(lmbda,G)
print("CM", cm)
c = c(K,M,cm) 
print("c",c)


class AnalyticalSol(UserExpression):
    def __init__(self,L,PL,c,Cm,M,Ku,G,alpha=1.0, t=1.0, **kwargs):
        super().__init__(**kwargs)
        self.L = L
        self.PL = PL
        self.c = c
        self.Cm = Cm
        self.t = t
        self.p0 = alpha*M*PL/(Ku + 4.*G/3.)
        self.uy_0_term = -PL/(Ku + 4.*G/3.)

    def eval(self,value, x):
        """ Evaluate the pressure expression """
        uy = 0.0
        p = 0.0
        for i in range(50):
            k = (2*i+1)*np.pi/(2*self.L)
            k2 = k*k
            expo_term = exp(-k2*self.c*self.t) 
            theta = k*(self.L+x[1])
            uy = uy +  (1./((2*i+1)**2))*expo_term*cos(theta)
            p = p + (1./(2*i + 1))*expo_term*sin(theta) 

        uy = self.Cm*self.p0*(-x[1] - (8.*L/(np.pi*np.pi))*uy) + self.uy_0_term*x[1]
        p = (4.*self.p0/np.pi)*p
        
        value[0] = 0        # ux
        value[1] = uy       # uy 
        value[2] = p        # pressure 
    
    def value_shape(self):
        return (3,)



mesh = RectangleMesh.create([Point(0.0,0.0), Point(1.0,10.0)],[nx,ny], CellType.Type.triangle) 

U = VectorElement("Lagrange", mesh.ufl_cell(),2, 2) # displacement 
P = FiniteElement("Lagrange", mesh.ufl_cell(),1) # pressure
W = FunctionSpace(mesh, U*P)

w = Function(W)
sol = AnalyticalSol(L,PL, c,cm,M,Ku,G, degree=2)
print("p0: ", sol.p0)

fileu = File("results/analytical/u_ana.pvd")
filep = File("results/analytical/p_ana.pvd")
    
day= 86400
time_steps = [0.1*day, day, 3*day, 5*day, 10*day]
#time_steps = [60,600,1800,3600]

for t_ in time_steps:
    sol.t = t_
    w.interpolate(sol)
    fileu << (w.split()[0],t_)
    filep << (w.split()[1],t_)
