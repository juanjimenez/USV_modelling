#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:07:56 2023

@author: juanjimenez
"""
import Bezier as Bz
import numpy as np
import math
import matplotlib.pyplot as pl

#creamos unos puntos de control de la curva de bezier
p = np.array([[1,3,1,3,1],[1,-4,0,4,-1]])
dt = 0.001 #paso del parametro p para construir las curvas
cbz,V,T = Bz.bezier(p,dt)

cbzp,Vp,Tp = Bz.bezier(p,100*dt)
pd, cdzd,Vd,Td = Bz.devbezier(p,100*dt)

maxx = max(p[0,:])
minx = min(p[0,:])
maxy = max(p[1,:])
miny = min(p[1,:])

pl.plot(cbz[0,:],cbz[1,:])
for i in range(cbzp.shape[1]):
    pl.arrow(cbzp[0,i],cbzp[1,i],cdzd[0,i],cdzd[1,i])
