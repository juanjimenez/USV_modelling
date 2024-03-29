#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:32:32 2024

@author: juan
"""
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.integrate import solve_ivp as sl
import control as ctrl

def USV(t,p,p0,A):
    dotp = A@p + p0 + np.array([0,2,0,0,0]) #pertubation added 
    return dotp
    
    
mu = np.array([[10.,0.],[0.,100.]])
mua = 100. #resistencia al giro 
mus = [mu,mua]
#_______calculo ganancias del controlador__________________________________
#control rumbo, control de la  fuerza posicion, contror error estimador
#calculamos primero las ganacias del estimador simplificado
A = np.array([[0,1],[0,-mu[0,0]]])
B = np.array([[0],[1]])
C = np.array([[1,0]])
Poles = [-10,-20] #perhaps improve it with LQR¿?
L = ctrl.acker(A.T,C.T,Poles).T
#y a continuación las del controlador, tambien simplificado (solo para p_x)
#incluimos acción integral para tratar de controlar las corrientes 
Pols = [-1,-2,-3]
Amp = np.array([[0, 1, 0],[0, -mu[0][0], 0], [-1, 0, 0]])
Bmp = np.array([[0],[1],[0]])
K = ctrl.acker(Amp,Bmp,Pols)

#
Atot = np.block([[A,  -B.dot(K[0,0:2].reshape(1,2)), -B.dot(K[0,2])],\
                 [L.dot(C), A-B.dot(K[0,0:2].reshape(1,2))-L.dot(C), -B.dot(K[0,2])],\
                 [-C, np.zeros(C.shape), np.zeros([C.shape[0],1])]]) 

# State vector
# p[0] ->> px
# p[1] ->> vx
# p[2] ->> ^px
# p[3] ->> ^vx
# p[4] ->> e 
pini = np.zeros(5)
pini[0] = 1.
pini[1] = 0.
pini[2] = 2.
pini[3] = 2.
pini[4] = 0
p0 = np.array([0,0,0,0,3])
tfin = 300. #2000.#tiempo final del experimento

#_______integrador LSODA (stiff)___________________________________________
USV(0,pini,p0,Atot)
sol = sl(USV, (0,tfin),pini,method='LSODA',args = (p0,Atot),max_step=0.05)

#dibujitos
pl.close('all')
pl.plot(sol.t,sol.y[0,:]) #posicion