# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 19:07:11 2022
Modelo del barco y algoritmo para estabilizarlo en un punto en presencia
de corrientes.
Es tan solo una prueba. Habría que refinar el mod en que se muestran los
resultados, etc.
Version modificada del original (estabilizacion) Para intentar emplear un
modelo que se coma el offset
@author: juan
"""

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.integrate import solve_ivp as sl
import vel_agua

campoux = vel_agua.campoux


#______Modelo del USV__________________________________________________________
def USV(t,p,E,mus,campvw,campo):
    """As a first approach we're going to include everything inside this 
    function
    t := time required by the integrations algorithm
    p := state variable vector
    E := [[0,-1],[1,0]] (pi/2 rotation matrix)
    mus :=(resistences) mu[0] body axis resistence to advance (Matrix). 
    mu[1] Resistence to turn around.
    camvm := function from vel_agua module to calculate the field of water 
    velocity
    campo := needed by vel_agua Model the water speed distribution
    '''
    State vector
    p[0] ->> px
    p[1] ->> py
    p[2] ->> vx
    p[3] ->> vy
    p[4] ->> R[0,0] (Rotation Matrix elements)
    p[5] ->> R[0,1]
    p[6] ->> R[1,0]
    p[7] ->> R[0,1]
    p[8] ->> w (angular speed)
    '''
    
    """

    ux = np.array([1,0])
    dotp = np.zeros(p.shape)
    mu = mus[0]
    mua = mus[1]
    wt = campvw(p[0:2],t ,campo)
    R = np.reshape(p[4:8],[2,2]) #USV heading 
    #torque to be applied open loop
    Tau = 15
    #FORCE TO BE APPLIED open loop
    F =  10 #Solo actuamos en la dir de avance
    
    #USV displacement
    dotp[0:2] = p[2:4] 
    dotp[2:4] = R@(ux*F) - R@(mu@(R.T@(p[2:4]- wt)))
    
    #Trusters force never used
    Td = (F + Tau)/2
    Ti = (F - Tau)/2
    Sw = E*p[8]
    dotR = Sw@R
    dotp[4:8] = np.reshape(dotR,[1,4])
    dotp[8] = Tau - mua*p[8] #desprecio cualquier par que pueda hacer la corriente
    
    
   
    return(dotp)
    

#_Inicializacion_______________________________________________________________    
'''
State vector
p[0] ->> px
p[1] ->> py
p[2] ->> vx
p[3] ->> vy
p[4] ->> R[0,0] (Rotation Matrix elements)
p[5] ->> R[0,1]
p[6] ->> R[1,0]
p[7] ->> R[0,1]
p[8] ->> w (angular speed)
'''
#thetai = np.arange(0,2*np.pi,np.pi/8) #initial USV heading
thetai = [np.pi/3]
pini = np.zeros(9)
r = 4 #distancia al origen a partir de la cual se prueba el algoritmo
ther = -np.pi/4 #Angulo posición respecto a eje x
v = 0.5 #v inicial del usv, (modulo)
#____________resolvemos para todos los angulos de thetai_______________________

for theta in thetai:
    pini[0] = r*np.cos(ther)
    pini[1] = r*np.sin(ther)
    pini[2] = v*np.cos(theta)
    pini[3] = v*np.sin(theta)
    pini[4] = np.cos(theta)
    pini[5] = -np.sin(theta)
    pini[6] = np.sin(theta)
    pini[7] = np.cos(theta)
    pini[8] = 0
    
    #Costantes del modelo absolutamente arbitrarias
    #matriz de rozamientos en ejes cuerpo. Incluyen todo Inercias, masas ..
    mu = np.array([[10.,0.],[0.,100.]])
    mua = 100. #resistencia al giro 
    mus = [mu,mua]
    #_______Campo de velocidades del agua _____________________________________
    E = np.array([[0,-1],[1,0]]) 
    campo = [np.pi/3,0,1,0.5,0] 
    tfin = 20. #2000.#tiempo final del experimento
    
    #_______integrador LSODA (stiff)___________________________________________
    USV(0,pini,E,mus,campoux,campo)
    sol = sl(USV,(0,tfin),pini,method='LSODA',args = (E,mus,campoux,campo),max_step=0.05)    
       
    
    #______Resultados (gráficos)_______________________________________________
    pl.close('all')
    pl.figure(1) #positions vs time
    pl.plot(sol.t,sol.y[0],'r')
    pl.plot(sol.t,sol.y[1],'k')
    pl.xlabel('t(s.)')
    pl.ylabel('m')
    pl.legend(['$p_1$','$p_2$'])
    pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    pl.savefig('./figures/pos_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
    pl.figure(2) #positions vs time
    pl.plot(sol.t,sol.y[2],'r')
    pl.plot(sol.t,sol.y[3],'k')
    pl.xlabel('t(s.)')
    pl.ylabel('m/s')
    pl.legend(['$\dot{p_1}$','$\dot{p_2}$'])
    pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    pl.savefig('./figures/vel_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
    pl.figure(3)
    pl.plot(sol.t,sol.y[8],'r')
    pl.xlabel('t(s.)')
    pl.ylabel('rad/s')
    pl.legend(['$\omega$'])
    pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    pl.savefig('./figures/omega_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
    
    
    vertices = np.array([[0, 0.25],[0,-0.25],[1,0],[0,0.25]])
    codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
    pl.figure(4) #trajectory and USV attitude
    
    pl.plot(sol.y[0],sol.y[1]) #drawing the trajetory
    
    #drawing the water speed field
    X = np.arange(-5,5,0.2)
    Y = np.arange(-5,5,0.2)
    vel_agua.pintacampoux(X,Y,campoux,campo,'g',tfin)
    L = sol.t.shape[0]
    
    #drawing the ship atitude
    step= 0.001
    
    ind = [int(step*((i-(L/2)**(1/3))**3+L/2)) for i in range(L)\
           if step*((i-(L/2)**(1/3))**3+L/2) < L]
        
    for i in ind:
        R = np.reshape(sol.y[4:8,i],[2,2])
        vertrot = np.array([np.dot(R,j) for j in vertices]) +\
        [sol.y[0,i],sol.y[1,i]]       
        pathi = Path(vertrot,codes)
        if i ==0:
            color = 'red'
        elif i == ind[-1]:
            color = 'green'
        else:
            color = 'blue'
        patchi = patches.PathPatch(pathi,facecolor = color)
        pl.gca().add_patch(patchi)
    
            
    
    pl.axis('equal')
    pl.xlabel('$p_1$(m.)')
    pl.ylabel('$p_2$(m.)')
    pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    #pl.savefig('./figures/disp_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
    
    # pl.figure(5)
    # pl.plot(sol.t,Td)
    # pl.plot(sol.t,Ti)
    # pl.legend(['$T_{S}$','$T_{P}$']) #P port babor izquierda, S starboard
    # pl.xlabel('t(s.)')
    # pl.ylabel('N.')
    # pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    # pl.savefig('./figures/thruster_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
