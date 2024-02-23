#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:55:47 2023
Ejemplo de camino paramétrico con curvas de bezier y representación de las
las superficies que conforman el camino
@author: juanjimenez
"""
import Bezier as Bz
import numpy as np
import math
import matplotlib.pyplot as pl

#creamos unos puntos de control de la curva de bezier
p = np.array([[1,3,1,3,1],[1,-4,0,4,-1]])
dt = 0.06 #paso del parametro p para construir las curvas
cbz,V,T = Bz.bezier(p,dt)
#buscamos los valores mayores y menores de los puntos de paso, para asegurar
#que la curva entra en la figura
maxx = max(p[0,:])
minx = min(p[0,:])
maxy = max(p[1,:])
miny = min(p[1,:])

x = np.linspace(minx,maxx,len(T))
y = np.linspace(miny,maxy,len(T))
xm,Tm = np.meshgrid(x,T)
xm,lmx = np.meshgrid(x,cbz[0,:])
ym,lmy = np.meshgrid(y,cbz[1,:])

ax = pl.figure().add_subplot(projection='3d')
ax.plot_wireframe(lmx,ym,Tm,linewidth=0.5,color='gray',label=r'$\phi_1(\xi)$')
ax.plot_wireframe(xm,lmy,Tm,linewidth=0.5,label=r'$\phi_2(\xi)$')
pl.plot(cbz[0,:],cbz[1,:],T,'r',linewidth=2.5,label=r'$^{aug}\mathcal{P}$')
pl.plot(cbz[0,:],cbz[1,:],np.zeros(T.shape),'b',linewidth=2.5,label=r'$\mathcal{P}$')
ax.legend()
ax.set_xlabel(r'$p_x$')
ax.set_ylabel(r'$p_y$')
ax.set_zlabel(r'$w$')