#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 17:17:48 2022

@author: juanjimenez
Marditos roedores. Estos son posibles campos para estudiar el efecto del cambio
de la velocidad del agua sobre la estabilidad del posicionamiento dinamico

"""

import numpy as np
import matplotlib.pyplot as pl

def campoux(p,t,campo):
   '''representa un campo uniforme en la direccion marcada por theta
   #el campo alcanza un maximo en la recta y = tan(theta)*x +b y cae a cero segu
   #nos alejamos de ella
   #p coordenadas en que se calcula el valor de la velocidad de arrastre
    el campo puede variar su orientacion con el tiempo si w no vale cero
   # A permite ajustar el valor de la velocidad máxima
   #como es una función gaussiana, S nos ajusta el  punto de inflexion
   #theta da la orientacion de la velocidad
   La variable campo es una lista.
   campo[0] => theta orientación (inicial) del campo
   campo[1] => w velocidad de rotación de la dirección del campo
   campo[3] => A velocidad máxima (modulo) del campo
   
   OJO HAY QUE REVISAR ESTO.
   '''
   theta = campo[0]
   w = campo[1]
   A = campo[2]
   S = campo[3]
   b = campo[4]
   #desrrotamos el punto. Solo necesito desrrotar la y
   pr = -p[0]*np.sin(theta+w*t) + p[1]*np.cos(theta+w*t)
   
   
   
   
   #el campo en el punto lo ligamos a     
   v=A*np.exp(-S*(pr+b)**2)

   wr = np.array([v*np.cos(theta+w*t),v*np.sin(theta+w*t)])
   return wr
   
def pintacampoux(X,Y,campoux,campo,col='b',t=0):
    '''para pintar el campo con flechillas 
    X e Y son dos vectores de puntos sobre los que hacer un grid de puntos en 
    los que calcular y pintar las flechillas del campo.
    campo es la misma lista que hemos definido en la funcion campoux
    '''
    
# =============================================================================
#     theta = campo[0]
#     w = campo[1]
#     A = campo[2]
#     S = campo[3]
#     b = campo[4]   
# =============================================================================
    Xm,Ym = np.meshgrid(X,Y)
    p =[Xm,Ym]
    Vm = campoux(p,0,campo)
    #Ymr = -Xm*np.sin(theta+w*t) + Ym*np.cos(theta+w*t)
    #Vm = A*np.exp(-S*(Ymr+b)**2)
    #Wmx = Vm*np.cos(theta+w*t)
    #Wmy = Vm*np.sin(theta+w*t)
    
    pl.quiver(Xm,Ym,Vm[0],Vm[1],color=col,width=1e-3)
            