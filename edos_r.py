import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Metodos_Chapra.edos import _fun_wrapped, _argumentos, RKF

def eulerIm(rhs, y0, t_eval, a):
    '''
    Parámetros
    ----------
    rhs : callable
        ecuación. Es necesario que como primer parámetro se pase la variable independiente (t, x, etcétera)
        
    y0 : lista o tupla
        condiciones iniciales.
        
    t_eval : ndarray
        malla de integración.
        
    a : float
        parte homogénea de la EDO


    Se considera el tipo de ecuación del estilo: f(t, y) = a*y(t) + q(t) 
    https://en.wikipedia.org/wiki/Backward_Euler_method
    
    
    Retorno
    -------
    X : ndarray
        solución a la ecuación o sistema de ecuaciones diferenciales.
    '''
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    for i, ti in enumerate(t_eval[1:]):
        yi = X[:, i]
        qt = _fun_wrapped(rhs, ti, yi) - a*yi
        X[:, i + 1] =  (yi + h*qt)/ (1 - a*h)
    
    return X


def heunIm(rhs, y0, t_eval):
    m, n, h = _argumentos(y0, t_eval)
    
    X = np.zeros((n, m)) ; X[:, 0] = y0
    
    #Tolerancia relativa
    atol = 1e-1
    
    #yi-1 = y0 ; yi = y1 y y2 = _predictor
    y1 = RKF(rhs, y0, t_eval)
    X[:, 1] = y1.T[1]
    
    def _predictor(f, t1, x0, x1, h):
        return x0 + _fun_wrapped(f, t1, x1)*2*h
    
    def _corrector(f, t1, t2, x1, x2, h):
        return x1 + (_fun_wrapped(f, t1, x1) + _fun_wrapped(f, t2, x2))/2 * h
    
    
    for i, ti in enumerate(t_eval[1:m - 1], start = 1):
        x1 = X[:, i] ; x0 = X[:, i - 1]
        t1 = ti
        t2 = t_eval[i + 1]
        x2 = _predictor(rhs, t1, x0, x1, h)
        imax = 100
        j = 0
        ea = 100
        
        while atol < ea and imax > j:
            #Corrector
            xj = _corrector(rhs, t1, t2, x1, x2, h)
            xjnew = xj
            xjold = x2
            ea = np.abs((xjnew - xjold)/xjnew)*100
            x2 = xjnew
            j += 1
            
        if j >= imax:
            raise ValueError('Iteraciones máximas alcanzada')
        
        X[:, i + 1] = xj
        
    return X




if __name__ == '__main__':
    DY_RIGIDA = lambda t, y : 4*np.exp(0.8*t) - 0.5*y  
    # DY_RIGIDA = lambda t, y : -1000*y + 3000 - 2000*np.exp(-t)
    IVC = 2
    T0 = 0 ; TF = 4
    METODOS = ['BDF', 'LSODA']
    n = 100 ; TMALLA = np.linspace(T0, TF, n)
    a = 0.5

    SOLivp = solve_ivp(DY_RIGIDA, [T0, TF], [IVC], METODOS[0], TMALLA)
    SOLeim = eulerIm(DY_RIGIDA, [IVC], TMALLA, a) #Acá gana el adaptativo por sobre RK4 aunque se densifique bastante la malla
    SOLhim = heunIm(DY_RIGIDA, [IVC], TMALLA)


    fig = plt.figure(figsize = (14, 7), dpi = 140)
    plt.scatter(TMALLA, SOLivp.y.T, color = 'green', label = 'IVP')
    plt.plot(TMALLA, SOLeim.T, color = 'blue', label = 'Euler Implicito')
    plt.plot(TMALLA, SOLhim.T, color = 'red', label = 'Heun sin autoinicio', linestyle = '--')
    plt.legend()
    