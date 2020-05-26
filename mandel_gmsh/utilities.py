import numpy as np
def Mc(phi,cf,cr=0,alpha=1.0):
    """Returns Biot Modulus (M) (Pa) [M L^-1 T^-2] in terms of compressibility

    :phi: porosity
    :cf: fluid compressibility Pa^-1 [M^-1 L T^2]
    :cr: rock grain compressibility  Pa^-1 [M^-1 L T^2]
    :alpha: Biot's coefficient
    :returns: M  

    """
    return 1./(phi*cf + (alpha-phi)*cr)


def Mb(phi=0.3,Kf=2.0e9,Ks=1.0e12,alpha=1.0):
    """Returns M, Biot modulus(Pa) [M L^-1 T^-2] in terms of compressibility

    :phi: porosity
    :Kf: bulk modulus of fluid [Pa] = 2.0e9 Pa
    :Ks: bulk modulus of rock grain/matrix [Pa] = 1.0e12 Pa
    :alpha: Biot's coefficient
    :returns: M (Pa)

    """
    return 1./(phi/Kf + (1.0-alpha)/Ks)


def Ku(lmbda, G, M, alpha=1.0):
    """Retruns Ku[Pa], undrained bulk modulus [M L^-1 T^-2]

    :lmbda: Lame parameter \lambda
    :G: Lame parameter \mu
    :M: Biot modulus
    :alpha: Biot coefficient
    :returns: Ku in smae unit of lmbda,G,M

    """
    return lmbda + 2*G/3. + alpha*alpha*M


def CM(lmbda, G):
    """Return vertical compressibility CM (Pa^-1)  [M^-1 L T^2]

    :lmbda: Lame parameter
    :G: Lame parameter
    :returns: CM

    """
    return 1./(lmbda + 2*G)



""" TODO: Fix the c method """
def c(K,M,CM, rho=1000, g=9.81, alpha=1.0 ):
    """Return consolidation coefficient c or cv [L^2 T^-1].

    :K: hydraulic conductivity tensor [m/s]
    :M: Compressibility [Pa]
    :CM: Vertical compressibility [Pa^-1]
    :rho: density of fluid = 1.0e3 [Kg/m^3]
    :g: gravity (9.81 [m/s^2]) 
    :alpha: Biot's coefficient
    :returns: c

    """
    tmp = rho*g*(1./M + alpha*alpha*CM)
    return (K/tmp)


def getB(M, Ku, alpha=1.0):
    """Return B,Skempton's pore pressure coefficient

    :M: Biot modulus
    :Ku: Undrained bulk modulus
    :alpha: Biot coefficient
    :returns: B

    """
    return alpha*M/(Ku)

def nu_u(nu,B,alpha=1.0):
    """Return undrainded Poisson's ratio

    :nu: Poisson's ratio
    :B: Skempton's pore pressure coefficient
    :alpha: Biot's coefficient
    :returns: nu_u

    """
    return (3.*nu + alpha*B*(1-2*nu))/(3.- alpha*B*(1-2*nu))



import sympy as smp
""" TODO: Fixthis. Returns wrong angles after 3 decimals"""
def getBeta(nu=0.25,nu_u=0.5):
    a = 0.0
    x = smp.Symbol('x')
    betas = []
    #coeff = nu/(nu_u - nu)
    coeff = (1.0 - nu)/(nu_u - nu)
    for j in range(20):
       for i in range(10):
           a = smp.nsolve(smp.sin(x) - coeff*x*smp.cos(x)/(a + i*(1-a)/9),a)
       betas.append(a)
       a += smp.pi.n() + 0.1
    betas[0] = 1.32760536836945
    betas = [2.0287578381103444, 4.913180439434637, 7.978665712413208, 11.085538406497022,14.207436725191188,
             17.33637792398336, 20.46916740274095, 23.604284772980407, 26.74091601478731,29.878586506107204,
             33.017001033357246, 36.15596641953655,39.295350981472986,42.43506188140989, 45.57503179559002,48.715210717557724,
             51.85556072915197, 54.99605255749639, 58.13666324489916,61.27737453356978,64.41817172183916]
    
    return betas





