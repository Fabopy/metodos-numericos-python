import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def error(xu, xi):
    return np.abs(xu - xi)/xu*100

class Cerrado:
    pass


def biseccion(func, xl, xu, es = None):
   if xl > xu:
       a = np.array([xu, xl])
       xu, xl = np.flip(a)
       
   if es is not None:
       es = es
   else:
       es = 1e-7
       
   ea = 100
    
   while es < ea:
       xr = (xl + xu)/2
       p = func(xl)*func(xr)
       if p > 0:
           xl = xr
           xnew = (xl + xu)/2
       else:
           xu = xr
           xnew = (xl + xu)/2
           
       ea = error(xnew, xr) 
       
   return xr, ea
               

def regula_falsi(func, xl, xu, es = None):
    if xl > xu:
        a = np.array([xu, xl])
        xu, xl = np.flip(a)
        
    if es is not None:
        es = es
    else:
        es = 1e-7
        
    ea = 100
    
    while es < ea:
        xr = xu - func(xu)*(xl - xu)/(func(xl) - func(xu))
        p = func(xl)*func(xr)
        if p > 0:
            xrold = xr
            xl = xrold
            xnew = xu - func(xu)*(xl - xu)/(func(xl) - func(xu))
            ea = error(xnew, xrold)
        else:
            xrold = xr
            xu = xrold
            xnew = xu - func(xu)*(xl - xu)/(func(xl) - func(xu))
            ea = error(xnew, xrold)
    
    return xr, ea
    

def regula_falsi_mod(func, xl, xu, es = None):
    pass
    




v1 = 850
th = np.deg2rad(15)
R = 287.05287 ; gamma = 1.4 ; T = 288.15
a = np.sqrt(gamma*R*T)
M1 = v1/a
F = lambda x : np.tan(th) - 2*(M1**2*np.sin(x)**2 - 1)/(np.tan(th)*(M1**2*(gamma + np.cos(2*x)) + 2))
xl = np.deg2rad(10) ; xu = np.deg2rad(40)
x_malla = np.linspace(xl, xu, 100)


fig = plt.figure(figsize = (14, 7), dpi = 140)
plt.grid(linestyle = '-.')
plt.title('Método gráfico')
plt.xlabel(r'$\theta_0$')
plt.ylabel(r'$F$')

plt.plot(x_malla, F(x_malla))

# plt.savefig('Metodos_Chapra/imagenes/RAICES/metodo_grafico.png', dpi = 300, bbox_inches="tight")

cr, _ = biseccion(F, xl = xl, xu = xu)
cr_rg, _ = regula_falsi(F, xl, xu)

#Verificación
if np.isclose(F(cr), 0):
    print(f'{cr} es la raíz de f')
