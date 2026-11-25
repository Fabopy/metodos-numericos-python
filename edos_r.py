import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Metodos_Chapra.edos import _fun_wrapped, _argumentos

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


if __name__ == '__main__':
    DY_RIGIDA = lambda t, y : -1000*y + 3000 - 2000*np.exp(-t)
    IVC = 0
    T0 = 0 ; TF = 4
    METODOS = ['BDF', 'LSODA']
    n = 100 ; TMALLA = np.linspace(T0, TF, n)
    a = -1000

    SOLivp = solve_ivp(DY_RIGIDA, [T0, TF], [IVC], METODOS[0], TMALLA)
    SOLeim = eulerIm(DY_RIGIDA, [IVC], TMALLA, a) #Acá gana el adaptativo por sobre RK4 aunque se densifique bastante la malla


    fig = plt.figure(figsize = (14, 7), dpi = 140)
    plt.plot(TMALLA, SOLivp.y.T, color = 'green', label = 'IVP')
    plt.scatter(TMALLA, SOLeim.T, color = 'blue', label = 'Euler Implicito')
    plt.legend()
    