import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt



def euler(rhs, y0, t_eval):
    '''
    rhs : función de lado derecho
    y0 : condición inicial
    t_span : intervalo de integración
    '''
    n = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    Y = np.zeros(n) ; Y[0] = y0
    for i in range(0, n - 1, 1):
        Y[i + 1] = Y[i] + rhs(t_eval[i], Y[i])*h
        
    return Y


def euler_sistema(rhs, y0, t_eval):
    '''
    rhs : lista de las ecuaciones diferenciales. Solamente resuelve un sistema de ecuaciones
    y0 : sequence de las condiciones iniciales
    t_eval : malla
    '''
    #Cantidad de puntos
    t_eval = np.asarray(t_eval)
    n = t_eval.size

    #Tamaño del paso
    h = t_eval[1] - t_eval[0] 
    
    # #Cantidad de funciones rhs
    # j = len(rhs)
    
    # Y = np.zeros((j, n))
    
    # for i in range(0, j, 1) :
    #     Y[i][0] = y0[i]
        
    rhs1 = rhs[0] ; rhs2 = rhs[1]
    
    Y1 = np.zeros(n) ; Y1[0] = y0[0]
    Y2 = np.zeros(n) ; Y2[0] = y0[1]
    
    for i in range(0, n - 1, 1):
        Y1[i + 1] = Y1[i] + rhs1(t_eval[i], Y1[i], Y2[i])*h
        Y2[i + 1] = Y2[i] + rhs2(t_eval[i], Y1[i], Y2[i])*h
        
    return Y1, Y2

def heun(rhs, y0, t_eval):
    '''
    rhs : función de lado derecho
    y0 : condición inicial
    t_span : intervalo de integración
    '''
    n = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    Y = np.zeros(n) ; Y[0] = y0
    
    for i in range(0, n - 1, 1):
        Y[i + 1] = Y[i] + rhs(t_eval[i], Y[i])*h
        y0p = Y[i + 1]
        ybar = (rhs(t_eval[i], Y[i]) + rhs(t_eval[i + 1], y0p))/2
        Y[i + 1] = Y[i] + ybar*h
    
    return Y
    
def RK2(rhs, y0, t_eval):
    '''
    rhs : función de lado derecho
    y0 : condición inicial
    t_span : intervalo de integración
    
    Método de Ralston (a_2 = 2/3)
    '''
    n = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    Y = np.zeros(n) ; Y[0] = y0
    for i in range(0, n - 1, 1):
        k1 = rhs(t_eval[i], Y[i])
        k2 = rhs(t_eval[i] + 3/4*h, Y[i] + 3/4*k1*h)
        Y[i + 1] = Y[i] + (1/3*k1 + 2/3*k2)*h
    
    return Y
    

def RK4(rhs, y0, t_eval):
    '''
    rhs : lista de las ecuaciones diferenciales. Solamente resuelve un sistema de ecuaciones
    y0 : condición inicial
    t_span : intervalo de integración
    '''
    n = t_eval.size
    h = t_eval[1] - t_eval[0] 
    
    rhs1 = rhs[0] ; rhs2 = rhs[1]
    
    Y1 = np.zeros(n) ; Y1[0] = y0[0]
    Y2 = np.zeros(n) ; Y2[0] = y0[1]
    
    for i in range(0, n - 1, 1):
        #k1
        k11 = rhs1(t_eval[i], Y1[i], Y2[i])
        k12 = rhs2(t_eval[i], Y1[i], Y2[i])
        
        #k2
        k21 = rhs1(t_eval[i] + 1/2*h, Y1[i] + 1/2*k11*h, Y2[i] + 1/2*k12*h)
        k22 = rhs2(t_eval[i] + 1/2*h, Y1[i] + 1/2*k11*h, Y2[i] + 1/2*k12*h)
        
        #k3
        k31 = rhs1(t_eval[i] + 1/2*h, Y1[i] + 1/2*k21*h, Y2[i] + 1/2*k22*h)
        k32 = rhs2(t_eval[i] + 1/2*h, Y1[i] + 1/2*k21*h, Y2[i] + 1/2*k22*h)
        
        #k4
        k41 = rhs1(t_eval[i] + h, Y1[i] + k31*h, Y2[i] + k32*h)
        k42 = rhs2(t_eval[i] + h, Y1[i] + k31*h, Y2[i] + k32*h)
        
        Y1[i+1] = Y1[i] + 1/6*(k11 + 2*k21 + 2*k31 + k41)*h
        Y2[i+1] = Y2[i] + 1/6*(k12 + 2*k22 + 2*k32 + k42)*h
        
    
    return Y1, Y2
        
 
if __name__ == '__main__':
    pass
    
    