"""
Created on Tue May 05 21:14:00 2026
Ajuste de un círculo a un conjunto de puntos
Para ajustar un círculo a un conjunto de puntos, podemos utilizar el método de mínimos cuadrados. El objetivo es encontrar el centro (h, k) y el radio r del círculo que minimicen la suma de las distancias al cuadrado desde cada punto al círculo.
El círculo se puede representar con la siguiente ecuación:
(x - h)^2 + (y - k)^2 = r^2 
Dado un conjunto de puntos (x_i, y_i), podemos definir la función de error que queremos minimizar:
E(h, k, r) = Σ((x_i - h)^2 + (y_i - k)^2 - r^2)^2
Podemos utilizar la función least_squares de scipy para encontrar los valores 
de h, k y r que minimicen esta función de error. Aquí hay un ejemplo de cómo hacerlo:
import numpy as np
from scipy.optimize import least_squares as ls
def error(params, x, y):
    h, k, r = params
    return ((x - h)**2 + (y - k)**2 - r**2)
# Ejemplo de uso
x = np.array([1, 2, 3, 4, 5])  
y = np.array([1, 2, 3, 4, 5])
initial_guess = [2, 2, 1]  # suposición inicial para el centro y el radio
result = ls(error, initial_guess, args=(x, y))
h, k, r = result.x
print(f"Centro: ({h}, {k}), Radio: {r}")
En este ejemplo, se define la función de error que calcula 
la diferencia entre la distancia al cuadrado desde cada punto al círculo y 
el radio al cuadrado. Luego, se utiliza least_squares para minimizar esta 
función de error y encontrar los parámetros del círculo 
que mejor se ajustan a los puntos dados.
JJC (todo este rollo lo ha escrito ChatGPT, no es mío).
"""

import numpy as np
from scipy.optimize import least_squares as ls
import matplotlib.pyplot as pl  

def error(params, x, y):
    h, k, r = params
    return ((x - h)**2 + (y - k)**2 - r**2)
#tomamos los puntos del terminas de ipython, para probar
x = sol.y[0]  # coordenadas x de los puntos
y = sol.y[1]  # coordenadas y de los puntos
#tomamos como suposicion inicial la media de los puntos como centro y la distancia media al centro como radio
h_initial = np.mean(x)
k_initial = np.mean(y)
r_initial = np.mean(np.sqrt((x - h_initial)**2 + (y - k_initial)**2))
initial_guess = [h_initial, k_initial, r_initial]
result = ls(error, initial_guess, args=(x, y))
h, k, r = result.x
print(f"Centro: ({h}, {k}), Radio: {r}")
#dibujamos el círculo ajustado junto con los puntos y la trayectoria del USV
theta = np.linspace(0, 2*np.pi, 100)
x_circle = h + r * np.cos(theta)
y_circle = k + r * np.sin(theta)
pl.figure()
pl.plot(x, y, 'ro', label='Puntos') 
pl.plot(x_circle, y_circle, 'b-', label='Círculo Ajustado')
pl.plot(sol.y[0], sol.y[1], 'g-', label='Trayectoria del USV')
pl.xlabel('x')
pl.ylabel('y')
pl.legend()
pl.title('Ajuste de Círculo a Puntos de Trayectoria del USV')
pl.axis('equal')
