import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def _fun_wrapped(func, t, y):
    return np.asarray(func(t, y))

def _argumentos(y0, t_eval):
    m = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    y0 = np.asarray(y0, dtype = float)
    n = len(y0) #Identifico la cantidad de ecuaciones que tiene mi sistema
    
    return m, n, h
    

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
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
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
    
    for i, ti in enumerate(t_eval[:-1]):
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
        for i, ti in enumerate(t_eval[:-1]):
            k1 = _fun_wrapped(rhs, ti, X[:, i])
            k2 = _fun_wrapped(rhs, ti + h, X[:, i] + k1*h)
            X[:, i + 1] = X[:, i] + (a2*k1 + a2*k2)*h
            
    if a2 == 1:
        for i, ti in enumerate(t_eval[:-1]):
            k1 = _fun_wrapped(rhs, ti, X[:, i])
            k2 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k1*h)
            X[:, i + 1] = X[:, i] + a2*k2*h
        
    if a2 == 2/3:
        for i, ti in enumerate(t_eval[:-1]):
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
    
    for i, ti in enumerate(t_eval[:-1]):
        k1 = _fun_wrapped(rhs, ti, X[:, i])
        k2 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k1*h)
        k3 = _fun_wrapped(rhs, ti + 1/2*h, X[:, i] + 1/2*k2*h)
        k4 = _fun_wrapped(rhs, ti + h, X[:, i] + k3*h)
        
        X[:, i + 1] = X[:, i] + 1/6*(k1 + 2*k2 + 2*k3 + k4)*h
        
    return X

def RKF(rhs, y0, t_eval):
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
    
    for i, ti in enumerate(t_eval[:-1]):
        k1 = _fun_wrapped(rhs, ti, X[:, i])
        k2 = _fun_wrapped(rhs, ti + 1/5*h, X[:, i] + 1/5*k1*h)
        k3 = _fun_wrapped(rhs, ti + 3/10*h, X[:, i] + 3/40*k1*h + 9/40*k2*h)
        k4 = _fun_wrapped(rhs, ti + 3/5*h, X[:, i] + h*(3/10*k1 - 9/10*k2 + 6/5*k3))
        k5 = _fun_wrapped(rhs, ti + h, X[:, i] + h*(-11/54*k1 + 5/2*k2 - 70/27*k3 + 35/27*k4))
        k6 = _fun_wrapped(rhs, ti + 7/8*h, X[:, i] + h*(1631/55296*k1 + 175/512*k2 + 575/13824*k3 + 44275/110592*k4 + 253/4096*k5))
        
        X[:, i + 1] = X[:, i] + (2825/27648*k1 + 18575/48384*k3 + 13525/55296*k4 + 277/14336*k5 + 1/4*k6)*h #acá en el chapra tiene un error, dice 18575/43384, pero en la figura 25.21 está bien escrito
        
    return X
    

def RK4A(rhs, y0, t_eval, atol = None):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación o sistema de ecuaciones diferenciales. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera) y luego el vector estado
        
    y0 : sequence
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.
        
    atol : Opcional
        Tolerancia absoluta. Como valor default es 1e-6

     Retorno
     -------
     X : ndarray
         solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    if atol is not None:
        atol = atol
    else:
        atol = 1e-6
        
    def MRK4(f, ti, yi, h):
        k1 = _fun_wrapped(f, ti, yi)
        k2 = _fun_wrapped(f, ti + 1/2*h, yi + 1/2*k1*h)
        k3 = _fun_wrapped(f, ti + 1/2*h, yi + 1/2*k2*h)
        k4 = _fun_wrapped(f, ti + h, yi + k3*h)
        
        return yi + 1/6*(k1 + 2*k2 + 2*k3 + k4)*h
        
    m, n, h = _argumentos(y0, t_eval)
 
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        yi = X[:, i]
        y1 = MRK4(rhs, ti, yi, h) #Predicción con un solo paso
        y2m = MRK4(rhs, ti, yi, h/2) #Predicción con dos medios pasos
        y2 = MRK4(rhs, ti + h/2, y2m, h/2) 
        d = y2 - y1
        y2 = y2 + d/15
        
        X[:, i + 1] = y2
        
    return X
            

if __name__ == '__main__':
    fig = plt.figure(figsize = (14, 7), dpi = 140)
    ax0, ax1 = fig.subplots(2, 1)
    
    def lotka_volterra(t, S):
        A = 1.5 ; B = 1 ; C = 3 ; D = 1
        x, y = S
        dx = A*x - B*x*y
        dy = -C*y + D*x*y
        
        dS = [dx, dy]
        return dS
    
    
    x0 = 10 ; y0 = 5
    t0 = 0 ; tf = 10
    n = 50 #Con esta cantidad de puntos se puede ver la ventaja de usar el RK4A
    t_malla = np.linspace(t0, tf, n)
    metodo = 'RK45'
    
    
    SOLivp = solve_ivp(lotka_volterra, t_span = [t0, tf], y0 = [x0, y0], t_eval = t_malla, method = metodo)
    
    SOLrk4 = RK4(lotka_volterra, y0 = [x0, y0], t_eval = t_malla)
    
    #Método adaptativo
    SOLrk4A = RK4A(lotka_volterra, y0 = [x0, y0], t_eval = t_malla)
    
    
    ax0.plot(t_malla, SOLivp.y[0], color = 'red', label = 'x(t)')
    ax0.plot(t_malla, SOLivp.y[1], color = 'blue', label = 'y(t)')

    ax0.plot(t_malla, SOLrk4A[0], color = 'red', label = 'x(t) con RK4 adaptativo', linestyle = (0, (1, 1)))
    ax0.plot(t_malla, SOLrk4A[1], color = 'blue', label = 'y(t) con RK4 adaptativo', linestyle = (0, (1, 1)))
    
    ax0.legend()
    
    
    ax1.plot(t_malla, SOLivp.y[0], color = 'red', label = 'x(t)')
    ax1.plot(t_malla, SOLivp.y[1], color = 'blue', label = 'y(t)')

    ax1.plot(t_malla, SOLrk4[0], color = 'red', label = 'x(t) con RK4', linestyle = '-.')
    ax1.plot(t_malla, SOLrk4[1], color = 'blue', label = 'y(t) con RK4', linestyle = '-.')
    
    
    ax1.legend()
    
    plt.tight_layout()
    
    
    

    
    