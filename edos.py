import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def _fun_wrapped(func, t, y):
    return np.asarray(func(t, y))

def euler(rhs, y0, t_eval):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación o sistema de ecuaciones diferenciales. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera) y luego el vector estado
        
    y0 : sequence
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.

    Retorno
    -------
    X : ndarray
        solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    m = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    y0 = np.asarray(y0, dtype = float)
    n = len(y0) #Identifico la cantidad de ecuaciones que tiene mi sistema
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(range(0, m - 1)):
        X[:, i + 1] = X[:, i] + h*_fun_wrapped(rhs, ti, X[:, i])
    
    return X


def heun(rhs, y0, t_eval):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación o sistema de ecuaciones diferenciales. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera) y luego el vector estado
        
    y0 : sequence
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.

    Retorno
    -------
    X : ndarray
        solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    m = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    y0 = np.asarray(y0, dtype = float)
    n = len(y0) 
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(range(0, m - 1)):
        X[:, i + 1] = X[:, i] + h*_fun_wrapped(rhs, ti, X[:, i])
        y0p = X[:, i + 1]
        ybar = (_fun_wrapped(rhs, ti, X[:, i]) + _fun_wrapped(rhs, t_eval[i + 1], y0p))/2
        X[:, i + 1] = X[:, i] + ybar*h
    
    return X
    
def RK2(rhs, y0, t_eval, a2 = None):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación o sistema de ecuaciones diferenciales. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera) y luego el vector estado
        
    y0 : sequence
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.
        
    a_2 : float or int para el caso de que a_2 = 1
     - a_2 = 1/2 Método de Heun con un solo corrector
     - a_2 = 1 Método del punto medio
     - a_2 = 2/3 Método de Ralston 
     
     Retorno
     -------
     X : ndarray
         solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    A2 = [1/2, 1, 2/3]
    
        
    if a2 not in A2:
        raise ValueError(f'a2 = {a2} no corresponde a ningún método')
    
    m = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    y0 = np.asarray(y0, dtype = float)
    n = len(y0)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    if a2 == 1/2:
        for i, ti in enumerate(range(0, m - 1)):
            k1 = _fun_wrapped(rhs, ti, X[:, i])
            k2 = _fun_wrapped(rhs, ti + h, X[:, i] + k1*h)
            X[:, i + 1] = X[:, i] + (a2*k1 + a2*k2)*h
            
    if a2 == 1:
        for i, ti in enumerate(range(0, m - 1)):
            k1 = _fun_wrapped(rhs, ti, X[:, i])
            k2 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k1*h)
            X[:, i + 1] = X[:, i] + a2*k2*h
        
    if a2 == 2/3:
        for i, ti in enumerate(range(0, m - 1)):
            k1 = _fun_wrapped(rhs, ti, X[:, i])
            k2 = _fun_wrapped(rhs, ti + 3/4*h, X[:, i] + 3/4*k1*h)
            X[:, i + 1] = X[:, i] + (1/3*k1 + a2*k2)*h
        
    return X
        

def RK4(rhs, y0, t_eval):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación o sistema de ecuaciones diferenciales. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera) y luego el vector estado
        
    y0 : sequence
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.
        
    a_2 : float or int para el caso de que a_2 = 1
     - a_2 = 1/2 Método de Heun con un solo corrector
     - a_2 = 1 Método del punto medio
     - a_2 = 2/3 Método de Ralston 
     
     Retorno
     -------
     X : ndarray
         solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    m = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    y0 = np.asarray(y0, dtype = float)
    n = len(y0)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(range(0, m - 1)):
        k1 = _fun_wrapped(rhs, ti, X[:, i])
        k2 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k1*h)
        k3 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k2*h)
        k4 = _fun_wrapped(rhs, ti + h, X[:, i] + k3*h)
        
        X[:, i + 1] = X[:, i] + 1/6*(k1 + 2*k2 + 2*k3 + k4)*h
        
    return X

 
if __name__ == '__main__':
    pass
    
    