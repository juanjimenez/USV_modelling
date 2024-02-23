# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 19:07:11 2022
Modelo del barco y algoritmo para estabilizarlo en un punto en presencia
de corrientes.
Es tan solo una prueba. Habría que refinar el mod en que se muestran los
resultados, etc.
@author: juan
"""

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
import vel_agua


pl.close('all')
#Costantes del modelo absolutamente arbitrarias
#matriz de rozamientos en ejes cuerpo. Incluyen todo Inercias, masas ..
mu = np.array([[10.,0.],[0.,100.]])
mua = 10. #resistencia al giro 

#ganancias del controlador
k1 = 3. #control de rumbo
k2 = 2.5 #control de la  fuerza (velocidad)
ki = 0.#01#01#control integral

InD = 0. #inicialización de la integral del error

tfin = 2000. #2000.#tiempo final del experimento
t = 0. #tiempo actual
dt = 0.01 #incremento de tiempos para integración

tp = tfin/150 #tiempo de pintar
P0i = np.array([[0.],[0.]]) #posicion de equilibrio
p = np.array([[1.],[-2.]]) #posición actual

#velocidad del agua, Camarón que se duerme...
#wt = np.array([[0.5],[-0.3]])
campo = [2*np.pi/3,0*np.pi/100,1,0.5]
wt = vel_agua.campoux(p, t, campo)
#R = np.eye(2)
#odio inicializar R asi pero no veo otra para tener cualquier orientación
theta = -np.pi/1
R = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
#Giro 90 natural 
E = np.array([[0.,-1.],[1.,0.]])

#velocidad angular del bicho
w = 0. 
#matriz skewsimetrica para construir \dot{R}
Sw = E*w

#fuerza de los thrusters
Td = 0.  #(Estribor)
Ti = 0.  #(babor)

#vector unitario en la direccion x (en ejes cuerpo es el heading del bicho)
ux = np.array([[1],[0]])

#pa pintar los resultados (triangulo que representa el barco)
vertices = np.array([[0, 0.25],[0,-0.25],[1,0],[0,0.25]])
codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]

#invento condiciones iniciales en las que supongo que me muevo en linea recta
#y el bicho ha alcanzado un estado estacionario 
F = 0 #10 #fuerza (fija la velocidad)
v = F*R.dot(np.linalg.inv(mu).dot(ux)) #velocidad actual
Tau = 0 #par

#pinto valores iniciales
# pl.figure(1)
# pl.plot(t,v[0],'.r')
# pl.plot(t,v[1],'.k')
# pl.plot(t,w,'.b')

# pl.figure(2)
# pl.plot(t,Tau,'.b')
# pl.plot(t,F,'.k')
# pl.plot(t,InD,'.r')


pl.figure(3)
pl.plot(p[0][0],p[1][0],'.')
#pintar un triangulo orientado en la posición del barco
vertrot = np.array([np.dot(R,j) for j in vertices]) +\
[p[0][0],p[1][0]]       
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'red')
pl.gca().add_patch(patchi)
theta = np.arange(0,2*np.pi+np.pi/50,np.pi/50)


# pl.figure(4)
# pl.plot(0,p[0][0],'.k')
# pl.plot(0,p[1][0],'.r')

# pl.figure(5)
# pl.plot(t,Ti,'r.')
# pl.plot(t,Td,'b.')


#bucle integración
#flanco= np.linalg.norm(P0-p) > 0.05
flag = 0

nD = np.linalg.norm(P0i-p)
P0 = P0i
ptotal =p
while t <= tfin:
    nD = np.linalg.norm(P0-p)
    Dd = (P0-p)/(nD+ (nD==0)) #rumbo al origen
    #para cambio de rumbo 
    st= Dd.T.dot(E.dot(R.dot(ux)))
    #para sentido de empuje
    ct = Dd.T.dot(R.dot(ux))
    if flag==1:
        if (abs(w)< 1e-6)&(np.linalg.norm(v)<1e-6):
            P0 = P0i + np.sign(ct)*R.dot(np.array([[nD],[0]]))
            flag = 0
            print(flag)
    else:
          if (abs(w)>=1e-3)&(np.linalg.norm(v)>=1e-6):
              P0 = P0i
              flag = 1
    
    # flancoold = flanco
    # flanco = int(nD>0.05)
    # if (flanco - flancoold) < 0:
    #     flag = 0
    #par a aplicar
    
    Tau = k1*st*np.sign(ct) 
    #fuerza a aplicar
    
    #InD += ct*flag*nD*dt #algo de integral para que se coma el offset de la wt
    F = k2*ct*nD**2 + ki*InD
    Td = (F + Tau)/2
    Ti = (F - Tau)/2
    p += v*dt
    ptotal = np.concatenate((ptotal,p),axis = 1)
    v += (R.dot(ux*(Td+Ti)) - R.dot(mu.dot(R.T.dot(v - wt))))*dt
    Sw = E*w 
    R += Sw.dot(R)*dt
    w += (Td-Ti - mua*w)*dt #desprecio cualquier par que pueda hacer la corriente
    
    
   
    #solo pinto de vez en cuando
    if t>=tp:
        #pl.figure(1)
        # pl.plot(t,v[0][0],'.r')
        # pl.plot(t,v[1][0],'.k')
        # pl.plot(t,w,'.b')
        # pl.figure(2)
        # pl.plot(t,Tau,'.b')
        # pl.plot(t,F,'.k')
        # pl.plot(t,InD,'.r')
        pl.figure(3)
        # pl.plot(p[0][0],p[1][0],'.')
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
        tp += tfin/200 #150
    
                 
    
    t = t + dt
    #pinto la última 
    if t>=tfin:
        #pl.figure(1)
        # pl.plot(t,v[0][0],'.r')
        # pl.plot(t,v[1][0],'.k')
        # pl.plot(t,w,'.b')
        # pl.figure(2)
        # pl.plot(t,Tau,'.b')
        # pl.plot(t,F,'.k')
        # pl.plot(t,InD,'.r')
        pl.figure(3)
        # pl.plot(p[0][0],p[1][0],'.')
        vertrot = np.array([np.dot(R,j) for j in vertices]) +\
        [p[0][0],p[1][0]]       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'green')
        pl.gca().add_patch(patchi)
        #pl.arrow(-p[1],p[0],v[1],v[0]) #NED
        
        # pl.figure(4)
        # pl.plot(t,p[0][0],'.k')
        # pl.plot(t,p[1][0],'.r')
        
        # pl.figure(5)
        # pl.plot(t,Ti,'r.')
        # pl.plot(t,Td,'b.')
    #calcula la velocidad del agua par la proxima iteracion
    wt = vel_agua.campoux(p, t, campo)
#leyendas y wt    
#pl.figure(1)
# pl.legend(['vx','vy','w'])
# pl.figure(2)
# pl.legend(['Tau','F'])
pl.figure(3)
pl.plot(ptotal[0,:],ptotal[1,:],'y')
pl.axis([-4,4,-4,4])
X = np.arange(-4,4,0.2)
Y = np.arange(-4,4,0.2)
vel_agua.pintacampoux(X,Y,campo,t)
# pl.figure(4)
# pl.legend(['x','y'])
# pl.figure(5)
# pl.legend(['Td','Ti'])
#pl.quiver(0,0,wt[0][0],wt[1][0])


#pl.arrow(-p[1],p[0],v[1],v[0]) 

