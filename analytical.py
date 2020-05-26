import numpy as np
from utilities import *
from fenics import *

MPa = 1e6
K = 1e-8  # [m/s]
phi = 0.375
#cf = 2.0e9*MPa  # [Pa]
Kf = 2.0e9*MPa  # [Pa]
cf = 1./Kf
alpha = 1.0  # Biot coefficient
PL = 1e5  # use load in [Pa] 
lmbda = 40*MPa  # Lame's coefficient [Pa]
G = 40*MPa  # \mu Lame's coefficient [Pa

M = Mb(phi,Kf)
print("Mb=", M)
print("Mc=",Mc(phi,cf))
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
        self.p0 = alpha*M*PL/(Ku + 4*G/3.)
        self.uy_0_term = -PL/(Ku + 4*G/3.)

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


L = 10.0
nx=4
ny=10
mesh = RectangleMesh.create([Point(0.0,0.0), Point(1.0,10.0)],[nx,ny], CellType.Type.triangle) 

U = VectorElement("Lagrange", mesh.ufl_cell(),2, 2) # displacement 
P = FiniteElement("Lagrange", mesh.ufl_cell(),1) # pressure
W = FunctionSpace(mesh, U*P)

w = Function(W)
sol = AnalyticalSol(L,PL, c,cm,M,Ku,G)

fileu = File("analytical/u_ana.pvd")
filep = File("analytical/p_ana.pvd")

t = 0
T = 3601
dt = 600

while ( t < T):

    w.interpolate(sol)
    fileu << (w.split()[0],t)
    filep << (w.split()[1],t)
    # update
    t = t+ dt;
    sol.t = t







#p0 = alpha*M*PL/(Ku + 4*G/3.)
#
#print("P0=", p0)
#
#class p(UserExpression):
#    def __init__(self,t,p0,c,L, **kwargs):
#        super().__init__(**kwargs)
#        self.t = t
#        self.p0 = p0
#        self.c = c
#        self.L = L
#        self.total = 0
#
#    def eval(self,value, x):
#        """ Evaluate the pressure expression """
#        self.total = 0.0
#
#        for i in range(20):
#            k = (2*i+1)*np.pi/(2*self.L)
#            k2 = k*k
#            self.total = self.total + (1./(2*i + 1))*exp(-k2*self.t)*sin(k*(L+x[1])) 
#            #self.total = self.total + (1./(2*i + 1))*exp(-k2*self.c*self.t)*sin(k*(L+x[1])) 
#
#        total= (4.*self.p0/np.pi)*self.total
#        #print("total=", self.total)
#        value[0] = total 
#    
#
#    def value_shape(self):
#        return ()
#
#P = FunctionSpace(mesh,"Lagrange",1)
##P = FunctionSpace(mesh,"DG",0)
#
#
#u = Function(P) 
#
#pressure = p(50.0,p0,c,L)
#u.interpolate(pressure)
#filep = File("analytical/p.pvd")
#filep << (u, 0)
#while ( t < T):
#    t = t+ dt;
#    pressure.t = t
#    u.interpolate(pressure)
#    filep << (u, t)
#




#
#
#code = '''
#class MyFunc : public Expression
#{
#public:
#
#  std::shared_ptr<MeshFunction<std::size_t> > cell_data;
#
#  MyFunc() : Expression()
#  {
#  }
#
#void eval(Array<double>& values, const Array<double>& x,
#          const ufc::cell& c) const
#  {
#       
#       total = 0.0
#       for(size_t i = 0; i < 20; i++ ){
#          double k1 = (2*i)*pi/(4*L);
#          double k2 = k1*k2;
#          total = total + 1./(2*i + 1)*exp(-k2*self.c*self.t)*sin(k*(L+x[1])) 
#       }
#        total= ((4./np.pi)*self.p0)*total
#        values[0] = total
#
#  }
#};'''
#
