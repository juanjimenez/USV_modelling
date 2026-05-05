#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:55:47 2023
Ejemplo de camino paramétrico con curvas de bezier y representación de las
las superficies que conforman el camino
@author: juanjimenez
"""
#import Bezier as Bz
import numpy as np
import math
import matplotlib.pyplot as pl

#creamos unos puntos de control de la curva de bezier

w = np.arange(0,2*np.pi+np.pi/20,np.pi/20)
cx = 2*np.cos(w)
cy = np.sin(2*w)

x = np.linspace(cx.min(),cx.max(),len(cx))
y = np.linspace(cy.min(),cy.max(),len(cy))


xm,Tm = np.meshgrid(x,w)
xm,lmx = np.meshgrid(x,cx)
ym,lmy = np.meshgrid(y,cy)

ax = pl.figure().add_subplot(projection='3d')
ax.plot_wireframe(lmx,ym,Tm,linewidth=0.5,color='gray',label=r'$\phi_1(\xi)$')
ax.plot_wireframe(xm,lmy,Tm,linewidth=0.5,label=r'$\phi_2(\xi)$')
pl.plot(cx,cy,w,'r',linewidth=2.5,label=r'$^{aug}\mathcal{P}$')
pl.plot(cx,cy,np.zeros(w.shape),'b',linewidth=2.5,label=r'$^{phys}\mathcal{P}$')
ax.legend()
ax.set_xlabel(r'$p_x$')
ax.set_ylabel(r'$p_y$')
ax.set_zlabel(r'$w$')