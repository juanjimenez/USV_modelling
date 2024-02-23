#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:13:03 2023

@author: juanjimenez
"""
import Bezier as Bz
import Bezier as Bz
import gvfbz as gvf
import vel_agua

import numpy as np
import math
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches

#creamos unos puntos de control de la curva de bezier###################
p = np.array([[10,30,10,30,10],[10,-40,0,40,-10]])
dt = 0.01 #paso del parametro p para construir las curvas
k = np.array([[500],[500]]) #ganancia parte perpendicular del campo

cbz,V,T = Bz.bezier(p,dt) #construimos la curva de bezier

#Campo inicial para s=0 s:= parámetro de la curva de Bezier############# 
steps = 30 #pasos del grid para calcular valores del campo para dibujar
x = np.linspace(0,31,steps)
y = np.linspace(-15,16,steps)
xm,ym = np.meshgrid(x,y) #puntos del plano sobre los que calcular vectores¿?
X= np.array([xm,ym])

s = 0
xi =gvf.bzxi(X, p, s, k) #campo
nxi = np.sqrt(xi[0]**2+xi[1]**2)
ps = Bz.bezierpnt(p, s)  #posición actual sobre la gráfica
ax =pl.figure(1).add_subplot()
qv = ax.quiver(xm,ym,xi[0]/nxi,xi[1]/nxi)
pl.plot(cbz[0],cbz[1])
pl.plot(ps[0],ps[1],'o')


#parametros del modelo del USV yc iniciales####################################
#Costantes del modelo del USV absolutamente arbitrarias
#matriz de rozamientos en ejes cuerpo. Incluyen todo Inercias, masas ..
mu = np.array([[10.,0.],[0.,100.]])
mua = 10. #resistencia al giro 

pusv = np.array([[0.],[0.]]) #posición actual

#Matriz de rotación R obtenida a partir del angulo de heading del USV
theta = -np.pi/1 #USV heading
R = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
#Giro 90 natural 
E = np.array([[0.,-1.],[1.,0.]])

#vector unitario en la direccion x (en ejes cuerpo es el heading del USV)
ux = np.array([[1],[0]])

#elijo unas c iniciales en las que supongo  estado estacionario el USv se 
#mueve con velocidad v = cte y en linea recta w = 0,

#descomentar si se quiere definir la fuerza como condicion inicial y optener la
#v estacionaria correspondiente 
# F = 0 #10 #fuerza (fija la velocidad)
# v = F*R.dot(np.linalg.inv(mu).dot(ux)) #velocidad actual
# Tau = 0 #par

#comentar si se quiere definir la fuerza como condicion inicial y obtener la 
#la v estacionaria correspondiente

vc = np.array([[0],[0]])#v en ejes cuerpo, dar valor solo a la primera coordenada

v = R@vc #v en ejes tierra
F = (mu@vc)[0,0] #fuerza de los thrusters
Tau = 0 #par generado por los thrusters
Td = (F + Tau)/2
Ti = (F - Tau)/2


#velocidad angular del bicho
w = 0. 
#matriz skewsimetrica para construir \dot{R}
Sw = E*w




#pa pintar los resultados (triangulo que representa el barco)
vertices = np.array([[0, 0.25],[0,-0.25],[1,0],[0,0.25]])
codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]


pl.plot(pusv[0][0],pusv[1][0],'.')
#pintar un triangulo orientado en la posición del barco
vertrot = np.array([np.dot(R,j) for j in vertices]) +\
[pusv[0][0],pusv[1][0]]       
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'blue')
pl.gca().add_patch(patchi)
theta = np.arange(0,2*np.pi+np.pi/50,np.pi/50)



#bucle de integración


tfin = 100. #2000.#tiempo final del experimento
t = 0. #tiempo actual
dt = 0.01 #incremento de tiempos para integración

tp = tfin/150 #tiempo de pintar

#velocidad del agua si hay corrientes
#velocidad del agua, Camarón que se duerme...
#wt = np.array([[0.5],[-0.3]])
campo = [2*np.pi/3,0*np.pi/100,1,0.5] #campo de velocidades
wt = vel_agua.campoux(pusv, t, campo) #velocidad del agua en la posición del barco

while t <= tfin:
   
    #par a aplicar
    
    Tau = 0 
    #fuerza a aplicar
    
    #InD += ct*flag*nD*dt #algo de integral para que se coma el offset de la wt
    F = 0
    Td = (F + Tau)/2
    Ti = (F - Tau)/2
    pusv += v*dt
    v += (R.dot(ux*(Td+Ti)) - R.dot(mu.dot(R.T.dot(v - wt))))*dt
    Sw = E*w 
    R += Sw.dot(R)*dt
    w += (Td-Ti - mua*w)*dt #desprecio cualquier par que pueda hacer la corriente
    
    
   
    #solo pinto de vez en cuando
    if t>=tp:
        # pl.figure(1)
        # pl.plot(t,v[0][0],'.r')
        # pl.plot(t,v[1][0],'.k')
        # pl.plot(t,w,'.b')
        # pl.figure(2)
        # pl.plot(t,Tau,'.b')
        # pl.plot(t,F,'.k')
        # pl.plot(t,InD,'.r')
        
        #pl.figure(3)
        pl.plot(p[0][0],p[1][0],'.')
        vertrot = np.array([np.dot(R,j) for j in vertices]) +\
        [p[0][0],p[1][0]]       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'blue')
        pl.gca().add_patch(patchi)
        #pl.arrow(-p[1],p[0],v[1],v[0]) #NED
        
        # pl.figure(4)
        # pl.plot(t,p[0][0],'.k')
        # pl.plot(t,p[1][0],'.r')
        
        # pl.figure(5)
        # pl.plot(t,Ti,'r.')
        # pl.plot(t,Td,'b.')
        
       
        tp += tfin/150
    t = t + dt
    #calcula la velocidad del agua par la proxima iteracion
    wt = vel_agua.campoux(pusv, t, campo)
  
vel_agua.pintacampoux(x,y,campo,t)
