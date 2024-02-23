#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:19:46 2023
Lyapunov theorem graphic
@author: juanjimenez
"""
import Bezier as bz
import numpy as np
import matplotlib.pyplot as pl
ax = pl.axes(projection='3d')
r = np.arange(0,3,0.25)
for i in r:
      
      p = i*np.array([[-1,-1,3,1,3,-1,-1],[1,-2,-4,0,4,3,1]])
      cbz,V,T = bz.bezier(p,0.001)
      if i  == r[4]:
          colorin ='blue'
          lin = 1
     
      else:
          colorin = 'black'
          lin = 0.5
          
      ax.plot3D(cbz[0],cbz[1],i**2*np.ones(np.shape(cbz[0])),linewidth=lin,color = colorin)
     
theta = np.arange(0,2*np.pi+0.01*2*np.pi,0.01*2*np.pi)
R = 2.5
x=R*np.cos(theta)
y=R*np.sin(theta)
z = np.zeros(np.shape(x))
#ax.plot3D(x,y,z)
alpha= 0.8
x=alpha*np.cos(theta)
y=alpha*np.sin(theta)
z = np.zeros(np.shape(x))
#ax.plot3D(x,y,z)

x=R*np.cos(theta)
y=R*np.sin(theta)
z = r[4]**2*np.ones(np.shape(x))
ax.plot3D(x,y,z,'r')
x=alpha*np.cos(theta)
y=alpha*np.sin(theta)
z = r[2]*np.ones(np.shape(x))
ax.plot3D(x,y,z,'g')
pl.xlim((-5,7))
pl.ylim((-5,7))
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$V(x)$')

pl.figure()
r = np.arange(0,3,0.25)
for i in r:
      
      p = i*np.array([[-1,-1,3,1,3,-1,-1],[1,-2,-4,0,4,3,1]])
      cbz,V,T = bz.bezier(p,0.001)
      if i  == r[4]:
          colorin ='blue'
          lin = 1
     
      else:
          colorin = 'white'
          lin = 0.5
          
      pl.plot(cbz[0],cbz[1],linewidth=lin,color = colorin)
R = 2.4
x=R*np.cos(theta)
y=R*np.sin(theta)
pl.plot(x,y,'r')
alpha= 0.75
x=alpha*np.cos(theta)
y=alpha*np.sin(theta)
pl.plot(x,y,'g')
pl.axis('equal')