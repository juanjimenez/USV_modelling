# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 19:07:11 2022
Modelito elemental de la boya. Lo que cuenta es el guiado
@author: abierto
"""

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
import gvf
pl.close('all')
#Costantes del modelo absolutamente arbitrarias
#matriz de rozamientos en ejes cuerpo. Incluyen todo Inercias, masas ..
mu = np.array([[10,0],[0,100]])
mua = 10. #resistencia al giro 

#parametros del controlador
ke = 0.5
kd = 40           

tfin = 60#tiempo final del experimento
t = 0. #tiempo actual
dt = 0.001 #incremento de tiempos para integración

#parametros de la elipse a seguir (gran guarrada)
a = 0#1 #semieje x
b = 5 #semieje y
alpha = 0 #orientación de la elipse
p0 = np.array([[5],[5]]) #centro de la elipse


tp = tfin/40 #tiempo de pintar
p = np.array([[1.],[0.]]) #posición actual Considero coordenadas NED

wt = np.array([[-0.3],[-0.3]]) #velocidad del agua, Camarón que se duerme...
R = np.eye(2)

#Giro 90 natural 
E = np.array([[0.,-1.],[1.,0.]])


w = 0. #velocidad angular del bicho
Sw = E*w

Td = 0. #fuerza de los thrusters
Ti = 0.

#Velocidad inicial, no parto del reposo porque gvf, necesita un p_dot distinto
#de cero
ux = np.array([[1],[0]])

#pa pintar los resultados
vertices = np.array([[0.05, 0],[-0.05,0],[0,0.2],[0.05,0]])
codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]

#bucle de integración
Fi = 10 #fuerza (fija la velocidad de consigna)
v = np.array([[0.],[Fi/mu[0,0]]]) #velocidad actual
Tau = 0 #par
e,n,H = gvf.circulo(p,p0,a)
 
pl.figure(1)
pl.plot(t,v[0],'.r')
pl.plot(t,v[1],'.k')
pl.plot(t,w,'.b')

pl.figure(2)
pl.plot(t,Tau,'.b')
pl.plot(t,Fi,'.k')
#pl.plot(t,Tau2,'.y')
#pl.plot(t,circulo.e,'*k')
pl.plot(t,e,'.r')

pl.figure(3)
pl.plot(-p[1][0],p[0][0],'.')
vertrot = np.array([np.dot(R,j) for j in vertices]) +\
[-p[1][0],p[0][0]]       
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'blue')
pl.gca().add_patch(patchi)
theta = np.arange(0,2*np.pi+np.pi/50,np.pi/50)
xc = a*np.cos(theta) - p0[1][0]
yc = a*np.sin(theta) + p0[0][0]
pl.plot(xc,yc)
pl.plot(-p0[1][0],p0[0][0],'^')
pl.figure(3)
pl.axis('equal')

#prueba con el control de héctor para circular
#circulo = gvfH.Path_gvf_circle(p0[0][0],p0[1][0], 2)
#orientacion inicial del vehiculo

while t <= tfin:
    #de momento solo controlamos rumbo y dejamos la v fija como partimos del
    #reposo deberíamos dejar que coja velocidad antes de lanzar el algoritmo..
    
    #e,n,H = gvf.elipse(a, b, alpha, p, p0)
    e,n,H = gvf.circulo(p,p0,a)
    m = R[:,1].reshape(2,1)
    Tau, dotXd, aux = gvf.gvf_control_boat(p, v, e, n, H, ke, kd, 1, m)
    
    #control de celeridad
    #F = np.log(e+1)*Fi
    F = (1.-np.exp(-e))*Fi
    #F=Fi
    Td = (F + Tau)/2
    Ti = (F - Tau)/2
    v += (R.dot(ux*(Td+Ti)) - R.dot(mu.dot(R.T.dot(v - wt))))*dt
    p += v*dt
    w += (Td-Ti - mua*w)*dt #desprecio cualquier par que pueda hacer la corriente
    Sw = E*w 
    R += R.dot(Sw)*dt
    
   
    #solo pinto de vez en cuando
    if t>=tp:
        pl.figure(1)
        pl.plot(t,v[0][0],'.r')
        pl.plot(t,v[1][0],'.k')
        pl.plot(t,w,'.b')
        pl.figure(2)
        pl.plot(t,Tau,'.b')
        pl.plot(t,F,'.k')
        #pl.plot(t,Tau2,'.y')
        #pl.plot(t,circulo.e,'*k')
        pl.plot(t,e,'.r')
        pl.figure(3)
        pl.plot(-p[1][0],p[0][0],'.')
        vertrot = np.array([np.dot(R,j) for j in vertices]) +\
        [-p[1][0],p[0][0]]       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'blue')
        pl.gca().add_patch(patchi)
        #pl.arrow(-p[1],p[0],v[1],v[0]) #NED
        pl.figure(4)
        pl.plot(t,aux,'r.')
        pl.plot(t,dotXd,'.b')
       
        tp += tfin/100
    t = t + dt
pl.figure(1)
pl.legend(['vx','vy','w'])
pl.figure(2)
pl.legend(['Tau','F','e'])

#pl.arrow(-p[1],p[0],v[1],v[0]) 

