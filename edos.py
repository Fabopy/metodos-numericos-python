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

def _phi(A, K):
    A = np.asarray(A, dtype = float)
    K = np.asarray(K, dtype = float)
    return np.dot(A, K)
    

def _euler(f, ti, yi, h):
    return yi + _fun_wrapped(f, ti, yi)*h

def _heun(f, ti, ti1, yi, h):
    yii = yi + f(ti, yi)*h
    ybar = (_fun_wrapped(f, ti, yi) + _fun_wrapped(f, ti1, yii))/2
    return yi + ybar*h

def _rk2(f, ti, yi, h, a2):
    a1 = 1 - a2
    p1 = 1/(2*a2)
    q11 = p1
    k1 = _fun_wrapped(f, ti, yi)
    k2 = _fun_wrapped(f, ti + p1*h, yi + q11*k1*h)
    phi = _phi([a1, a2], [k1, k2])
    return yi + phi*h

def _rk4(f, ti, yi, h):
    a1 = 1/6 ; a2 = 1/3 ; a3 = 1/3 ; a4 = 1/6
    k1 = _fun_wrapped(f, ti, yi)
    k2 = _fun_wrapped(f, ti + 1/2*h, yi + 1/2*k1*h)
    k3 = _fun_wrapped(f, ti + 1/2*h, yi + 1/2*k2*h)
    k4 = _fun_wrapped(f, ti + h, yi + k3*h)
    phi = _phi([a1, a2, a3, a4], [k1, k2, k3, k4])
    
    return yi + phi*h

def _rkf(f, ti, yi, h):
    a1 = 2825/27648 ; a2 = 0 ; a3 = 18575/48384 ; a4 = 13525/55296 ; a5 = 277/14336 ; a6 = 1/4
    k1 = _fun_wrapped(f, ti, yi)
    k2 = _fun_wrapped(f, ti + 1/5*h, yi + 1/5*k1*h)
    k3 = _fun_wrapped(f, ti + 3/10*h, yi + 3/40*k1*h + 9/40*k2*h)
    k4 = _fun_wrapped(f, ti + 3/5*h, yi + h*(3/10*k1 - 9/10*k2 + 6/5*k3))
    k5 = _fun_wrapped(f, ti + h, yi + h*(-11/54*k1 + 5/2*k2 - 70/27*k3 + 35/27*k4))
    k6 = _fun_wrapped(f, ti + 7/8*h, yi + h*(1631/55296*k1 + 175/512*k2 + 575/13824*k3 + 44275/110592*k4 + 253/4096*k5))
    phi = _phi([a1, a2, a3, a4, a5, a6], [k1, k2, k3, k4, k5, k6])
    
    return yi + phi*h

def _rka(f, ti, yi, h):
    pass





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
        yi = X[:, i]
        X[:, i + 1] = _euler(rhs, ti, yi, h)
    
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
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        tii = t_eval[i + 1]
        yi = X[:, i]
        X[:, i + 1] = _heun(rhs, ti, tii, yi, h)
    
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
    
    if a2 is None:
        a2 = 2/3
    
    if a2 not in A2:
        raise ValueError(f'a2 = {a2} no corresponde a ningún método')
        
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        yi = X[:, i]
        X[:, i + 1] = _rk2(rhs, ti, yi, h, a2)
        
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
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        yi = X[:, i]
        X[:, i + 1] = _rk4(rhs, ti, yi, h)
        
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
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        yi = X[:, i]
        X[:, i + 1] = _rkf(rhs, ti, yi, h)
        
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
        
    m, n, h = _argumentos(y0, t_eval)
 
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[:-1]):
        yi = X[:, i]
        y1 = _rk4(rhs, ti, yi, h) #Predicción con un solo paso
        y2m = _rk4(rhs, ti, yi, h/2) #Predicción con dos medios pasos 
        y2 = _rk4(rhs, ti + h/2, y2m, h/2) #Predicción con dos medios pasos 
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
    
    
    

    
    