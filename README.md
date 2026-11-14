# Métodos Numéricos para Ingenieros en Python

Se utiliza el libro de Steven Chapra y Raymond Canale. 

Libro

[Numerical Methods for Engineers](https://www.amazon.com/Package-Numerical-Methods-Engineers-Semester/dp/1259289168/ref=sr_1_23?s=books&sr=1-23)

Los temas son:

1. [Raíces de ecuaciones](#raíces-de-ecuaciones)
2. [Ecuaciones Algebraicas Lineales](#ecuaciones-algebraicas-lineales)
3. [Optimización](#optimización)
4. [Ajuste de Curvas](#ajuste-de-curvas)
5. [Diferenciación e Integración Numérica](#diferenciación-e-integración-numérica)
6. [Ecuaciones Diferenciales Ordinarias](#ecuaciones-diferenciales-ordinarias)
7. [Ecuaciones Diferenciales Parciales](#ecuaciones-diferenciales-parciales)

## Raíces de ecuaciones

## Ecuaciones Algebraicas Lineales

## Optimización

## Ajuste de Curvas

## Diferenciación e Integración Numérica

## Ecuaciones Diferenciales Ordinarias

Este capítulo se dedica a la solución de ecuaciones diferenciales ordinarias de la forma

$$
\frac{dy}{dx} = f(x, y)
$$

De la definición de diferencias finitas hacia adelante, podemos despejar $y_{i+1}$

$$
y_{i+1} = y_i + \frac{dy}{dx}(x_{i+1} - x_i)
$$

donde generalmente a $dy/dx$ se le denota con $\phi$ y $(x_{i+1} - x_i) = h$ es el paso. Existen varias formas de calcular $y_{i+1}$ pero las más conocidas son:

- Método de Euler
- Método de Heun
- Runge-Kutta de 2do orden
- Runge-Kutta de 4to orden


### Método de Euler
La primera derivada ofrece una estimación directa de la pendiente en $x_i$:

$$
\phi = f(x_i, y_i)
$$


donde $f(x_i, y_i)$ es la ecuación diferencial evaluada en $x_i$ y $y_i$. La estimación se sustituye en la ecuación:

$$
y_{i+1} = y_i + f(x_i, y_i)(x_{i+1} - x_i)
$$

Esta fórmula se conoce como **método de Euler**. Veamos un ejemplo:

$$
\frac{dv}{dt} = g - \frac{c}{m}v
$$

La solución analítica de esta ecuación diferencial es:

$$
v(t) = \frac{gm}{c}\left(1 - e^{-c/m \ t}\right)
$$


![metodoEuler](C:/Users/fabri/OneDrive/Documents/PythonSpyder/foto_euler.pdf)





## Ecuaciones Diferenciales Parciales