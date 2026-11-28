import numpy as np
import matplotlib.pyplot as plt
from scipy.differentiate import derivative

def Richardson(D1, D2):
    raise NotImplementedError
    # return 4/3*D2 - 1/3*D1

def d_func4(f, x, i, h):
    return (-f(x[i + 2]) + 8*f(x[i + 1]) - 8*f(x[i - 1]) + f(x[i - 2]))/(12*h)

def d_func2(f, x, i, h):
    return (f(x[i + 1]) - f(x[i - 1]))/(2*h)


def d_funcInterpolacion(f, x, I):
    '''
    Se sigue la siguiente fórmula:
        a = x_{i - 1}
        b = x_i
        c = x_{i + 1}
        
        f'(x) = f(a) * (2*x - b - c) / [(a - b)*(a - c)] + f(b) * (2*x - a - c) / [(b - a)*(b - c)] + f(c)*(2*x - a - b)/[(c - a)*(c - b)]
    '''
    a, b, c = I
    A = (2*x - b - c)/((a - b)*(a - c))
    B = (2*x - a - c)/((b - a)*(b - c))
    C = (2*x - a - b)/((c - a)*(c - b))
    
    return f(a) * A + f(b) * B + f(c) * C




#Diferencias finitas centradas
def derO4(f, x):
    n = x.size
    f1 = np.array([])
    h1 = x[1] - x[0]
    for i in range(2, n - 2, 1):
        D = d_func4(f, x, i, h1)
        f1 = np.append(f1, D)
        
    malla = x[2:-2]
    
    return malla, f1

def derO2(f, x):
    n = x.size
    f1 = np.array([])
    h1 = x[1] - x[0]
    for i in range(1, n - 1, 1):
        D = d_func2(f, x, i, h1)
        f1 = np.append(f1, D)
        
    malla = x[1:-1]
    
    return malla, f1
    
def derInterpolacion(f, x):
    n = x.size
    f1 = np.array([])
    for i in range(1, n - 1, 1):
        I = (x[i - 1], x[i], x[i + 1])
        D = d_funcInterpolacion(f, x[i], I)
        f1 = np.append(f1, D)
        
    malla = x[1:-1]
    
    return malla, f1


if __name__ == '__main__':
    
    L = 4
    t = np.array([0, 0.1, 0.2, 0.3, 0.5, 0.7])
    i = np.array([0, 0.16, 0.32, 0.56, 0.84, 2])
    
    fig = plt.figure(figsize = (14, 7), dpi = 140)
    plt.grid(linestyle = '--')
    plt.xlabel('t')
    plt.ylabel(r'$V_L$')
    
    plt.plot()

