import numpy as np
from scipy.optimize import OptimizeResult

class RootResult(OptimizeResult):
    xr : float
    ea : float
    i : float
    
def error(xu, xi):
    return np.abs(xu - xi)/xu*100

class Cerrado:
    '''
    Clase para encontrar raíces mediante métodos cerrados. 
    
    Parámetros:
        func : función a encontrar la raíz
        xl : valor ínfimo del intervalo considerado
        xu : valor superior del intervalo considerado
    '''
    def __init__(self, func, xl, xu, es = None):
        self.f = func
        self.xl = xl
        self.xu = xu
        if es is not None:
            self.es = es
        else:
            self.es = 1e-7
            
                
    def biseccion(self):
        metodo = 'Bisección'
        if self.xl > self.xu:
            a = np.array([self.xu, self.xl])
            self.xu, self.xl = np.flip(a)
            

        ea = 100
        i = 0
        l = self.xl
        u = self.xu
        
        if self.f(l)*self.f(u) > 0:
            raise ValueError(f'El intervalo [{l}, {u}] no encierra a la raíz')
        
        while self.es < ea:
            xr = (l + u)/2
            p = self.f(l)*self.f(xr)
            if p > 0:
                l = xr
                xnew = (l + u)/2
            else:
                u = xr
                xnew = (l + u)/2
                
            ea = error(xnew, xr) 
            i += 1
        
        return RootResult(xr = xr, 
                          ea = ea,
                          i = i,
                          method = metodo)
    
    def regula_falsi(self):
        metodo = 'Falsa posición'
        if self.xl > self.xu:
            a = np.array([self.xu, self.xl])
            self.xu, self.xl = np.flip(a)
            
        g = lambda l, u : u - self.f(u)*(l - u)/(self.f(l) - self.f(u))
        
        ea = 100
        i = 0
        l = self.xl
        u = self.xu
        
        if self.f(l)*self.f(u) > 0:
            raise ValueError(f'El intervalo [{l}, {u}] no encierra a la raíz')
        
        while self.es < ea:
            xr = g(l, u)
            p = self.f(l)*self.f(xr)
            if p > 0:
                xrold = xr
                l = xrold
                xnew = g(l, u)
            else:
                xrold = xr
                u = xrold
                xnew = g(l, u)
                
            ea = error(xnew, xrold)
            i += 1
        
        return RootResult(xr = xr, 
                          ea = ea,
                          i = i,
                          method = metodo)
    
    def regula_falsi_mod(self):
        metodo = 'Falsa posición modificada'
        if self.xl > self.xu:
            a = np.array([self.xu, self.xl])
            self.xu, self.xl = np.flip(a)
        
        ea = 100
        i = 0
        l = self.xl
        u = self.xu
        fl = self.f(l)
        fu = self.f(u)
        il = 0
        iu = 0
        
        if self.f(l)*self.f(u) > 0:
            raise ValueError(f'El intervalo [{l}, {u}] no encierra a la raíz')
        
        while self.es < ea:
            xr = u - fu*(l - u)/(fl - fu)
            xrold = xr
            fr = self.f(xr)
            p = fl*fr
            if p < 0:
                u = xr
                fu = self.f(u)
                il = il + 1
                iu = 0
                xnew = u - fu*(l - u)/(fl - fu)
                ea = error(xnew, xrold)
                if il >= 2:
                    fl = fl/2
            else:
                l = xr
                fl = self.f(l)
                iu = iu + 1
                il = 0
                xnew = u - fu*(l - u)/(fl - fu)
                ea = error(xnew, xrold)
                if iu >= 2:
                    fu = fu/2
            i += 1
                    
        return RootResult(xr = xr, 
                          ea = ea,
                          i = i,
                          method = metodo)
            
                

        
class Abierto:
    def __init__(self, func, x0, es = None, imax = None):
        self.f = func
        self.x0 = x0
        if es is not None:
            self.es = es
        else:
            self.es = 1e-7
            
        if imax is not None:
            self.imax = imax
        else:
            self.imax = 500
            
            
    def PuntoFijo(self, g):
        '''
        Considerar que la forma de construir g(x) puede no converger y otras sí
        '''
        metodo = 'Punto Fijo'
        ea = 100
        i = 0
        xrold = self.x0
        
        while self.es < ea and self.imax > i:
            xr = g(xrold)
            xnew = xr
            ea = error(xnew, xrold)
            xrold = xnew
            i += 1
            
        if i >= self.imax:
            mensaje = 'Iteración máxima alcanzada'
            return RootResult(message = mensaje, 
                              xr = xr,
                              ea = ea,
                              i = i,
                              method = metodo)
        
        return RootResult(xr = xr,
                          ea = ea,
                          i = i,
                          method = metodo)
    
    def Newton_Raphson(self, dfunc):
        metodo = 'Newton Raphson'
        ea = 100
        i = 0
        xrold = self.x0
        tol = 1e-20
        
        while self.es < ea and self.imax > i:
            if np.abs(dfunc(xrold)) < tol:
                raise ZeroDivisionError(f'División por cero debido a que la derivada es cada vez más chica en magnitud: {dfunc(xrold)}')
            xr = xrold - self.f(xrold)/dfunc(xrold)
            xnew = xr
            ea = error(xnew, xrold)
            xrold = xnew
            i += 1
            
        return RootResult(xr = xr, 
                          ea = ea,
                          i = i,
                          method = metodo)
            
        
    def Secante(self, xb):
        metodo = 'Secante'
        ea = 100
        i = 0
        xrold = self.x0
        
        while self.es < ea and self.imax > i:
            xr = xrold - self.f(xrold)*(xb - xrold)/(self.f(xb) - self.f(xrold))
            xnew = xr
            ea = error(xnew, xrold)
            xb = xrold
            xrold = xnew
            i += 1
            
        return RootResult(xr = xr, 
                          ea = ea,
                          i = i,
                          method = metodo)
        
    def Secante_Modificada(self):
        pass
        
            





    
    



