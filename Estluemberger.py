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
import control as ctrl
import vel_agua

campoux = vel_agua.campoux


#______Modelo del USV__________________________________________________________
def USV(t,p,p0,mus,Ks,C,E,campvw,campo):
    """As a first approach we're going to include everything inside this 
    function
    t := time required by the integrations algorithm
    p := state variable vector
    p0 := target stabilised final position
    E := [[0,-1],[1,0]] (pi/2 rotation matrix)
    mus :=(resistences) mu[0] body axis resistence to advance (Matrix). 
    mu[1] Resistence to turn around.
    Ks := Controller constants k[0] Heading control k[1] position crtl
    camvm := function from vel_agua module to calculate the field of water 
    velocity
    campo := needed by vel_agua Model the water speed distribution
    T :truster initial values
    """

    ux = np.array([1,0])
    dotp = np.zeros(p.shape)
    mu = mus[0]
    mua = mus[1]
    kr= Ks[0]
    k = Ks[1]
    L = Ks[2]
    wt = campvw(p[0:2],t ,campo)
    R = np.reshape(p[8:12],[2,2]) #USV heading 
    ep =p0-p[0:2] #position error p0 is the desired equilibrium point to down the probe 
    epph = p[0:2] - p[4:6] #error of estimation position
    nD = np.linalg.norm(ep)
    Dd = ep/(nD+ (nD==0)) #heading towards the origin (equilibrium point)
    #
    st= Dd.T@E@R@ux
    #
    ct = Dd.T@R@ux
    #ev = p[2:4]-p[6:8] #velocity error
    #torque to be applied
    Tau = kr*st*(np.sign(ct))#+1*(ct==0))#la condición es para evitar el punto de eq inestable
    #Tau = kr*st/ct
    
    
    
    
    F =  k*(R.T@ep)[0] -(R@mu@R.T@p[13:15])[0]#Solo actuamos en la dir de avance
    
    #USV displacement
    dotp[0:2] = p[2:4] 
    dotp[2:4] = R@(ux*F) - R@(mu@(R.T@(p[2:4]- wt)))
    #dotp[2:4] = R.dot(ux*F) - R.dot(mu.dot(R.T.dot(p[2:4])))
    #USV rotation
    Td = (F + Tau)/2
    Ti = (F - Tau)/2
    Sw = E*p[12]
    dotR = Sw@R
    dotp[8:12] = np.reshape(dotR,[1,4])
    dotp[12] = Td-Ti - mua*p[12] #desprecio cualquier par que pueda hacer la corriente
    
    #USV estimator y estimador de corrientes
    
    
    dotp[4:6] = p[6:8] - L[2:4,:]@epph
    dotp[6:8] = R@(ux*F) - R@mu@R.T@(p[6:8]-p[13:15]) - L[0:2,:]@epph
    dotp[13:15] = L[4:6,:]@epph
    
    return(dotp)
    

#_Inicializacion_______________________________________________________________    
'''
State vector
p[0] ->> px
p[1] ->> py
p[2] ->> vx
p[3] ->> vy
p[4] ->> ^px (^. Estimator variables)
p[5] ->> ^py
p[6] ->> ^vx
p[7] ->> ^vy
p[8] ->> R[0,0] (Rotation Matrix elements)
p[9] ->> R[0,1]
p[10] ->> R[1,0]
p[11] ->> R[0,1]
p[12] ->> w (angular speed)
p[13] ->> ^vwx (estimated water speed)
p[14] ->> ^vwy (Estimated water speed)
'''
#thetai = np.arange(0,2*np.pi,np.pi/8) #initial USV heading
thetai = [np.pi/3]
pini = np.zeros(14)
r = 4 #distancia al origen a partir de la cual se prueba el algoritmo
v = 0.5 #v inicial del usv, (modulo)
#____________resolvemos para todos los angulos de thetai_______________________ 
for theta in thetai:
    pini[0] = r*np.cos(theta)
    pini[1] = r*np.sin(theta)
    pini[2] = -v*np.cos(theta)
    pini[3] = -v*np.sin(theta)
    pini[4] = r*np.cos(theta)
    pini[5] = r*np.sin(theta)
    pini[6] = -v*np.cos(theta)
    pini[7] = -v*np.sin(theta)
    pini[8] = np.cos(theta)
    pini[9] =  -np.sin(theta)
    pini[10] = np.sin(theta)
    pini[11] = np.cos(theta)
    pini[12] = 0
    pini[13] = 0
    pinni[14] = 0
    p0 = np.array([0,0]) #desired stabilization position
    
    #Costantes del modelo absolutamente arbitrarias
    #matriz de rozamientos en ejes cuerpo. Incluyen todo Inercias, masas ..
    mu = np.array([[10.,0.],[0.,100.]])
    mua = 100. #resistencia al giro 
    mus = [mu,mua]
    
    #_______calculo ganancias del controlador__________________________________
    #control rumbo, control de la  fuerza posicion, contror error estimador
    #calculamos primero las ganacias del estimador simplificado
    A = np.array([[0,1],[0,-mu[0,0]]])
    B = np.array([[0],[1]])
    C = np.array([[1,0]])
    Poles = [-3,-3] #perhaps improve it with LQR¿?
    L = ctrl.acker(A.T,C.T,Poles).T
    #y a continuación las del controlador, tambien simplificado (solo para p_x)
    #incluimos acción integral para tratar de controlar las corrientes 
    Pols = [-2,-2,-1]
    Amp = np.array([[0, 1, 0],[0, -mu[0][0], 0], [-1, 0, 0]])
    Bmp = np.array([[0],[1],[0]])
    K = 10 #de momento la ganancia proporcional al error en x la metemos a capon
    #la  ganancia del control de rotación está metida a pelo
    
    Ks = [50.,K,L] #
    
    E = np.array([[0,-1],[1,0]])
    #campo = [3*np.pi/4,0*np.pi/100,1,0.5]
    campo = [np.pi/3,0*np.pi/100,1,0.5]
    tfin = 20. #2000.#tiempo final del experimento
    
    #_______integrador LSODA (stiff)___________________________________________
    USV(0,pini,p0,mus,Ks,C,E,campoux,campo)
    sol = sl(USV, (0,tfin),pini,method='LSODA',args = (p0,mus,Ks,C,E,campoux,campo),max_step=0.05)    
       
    
    #______Resultados (gráficos)_______________________________________________
    pl.close('all')
    pl.figure(1) #positions vs time
    pl.plot(sol.t,sol.y[0],'r')
    pl.plot(sol.t,sol.y[1],'k')
    pl.plot(sol.t,sol.y[4],'g')
    pl.plot(sol.t,sol.y[5],'b')
    
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
    # # pl.plot(t,w,'.b')
    # pl.figure(2) #velocities vs time
    # pl.plot(sol.t,sol.y[2],'.r')pl.figure(1) #positions vs time
    pl.figure(3)
    pl.plot(sol.t,sol.y[12],'r')
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
    vel_agua.pintacampoux(X,Y,campo,tfin)
    L = sol.t.shape[0]
    
    #drawing the ship atitude
    step= 0.001
    
    ind = [int(step*((i-(L/2)**(1/3))**3+L/2)) for i in range(L)\
           if step*((i-(L/2)**(1/3))**3+L/2) < L]
        
    for i in ind:
        R = np.reshape(sol.y[8:12,i],[2,2])
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
    pl.savefig('./figures/disp_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
    
    ux = np.array([1,0])
    dotp = np.zeros(13)
    mu = mus[0]
    mua = mus[1]
    k1 = Ks[0]
    k2 = Ks[1]
    k3 = Ks[2]
   
    Td = np.zeros(L)
    Ti = np.zeros(L)
    for j in range(0,L):
        sly =sol.y[:,j]
        R = np.reshape(sly[8:12],[2,2])
        ep = p0-sly[0:2] #position error
    
        nD = np.linalg.norm(ep)
        Dd = ep/(nD+ (nD==0)) #heading towards the origin
    
        st= Dd.T.dot(E.dot(R.dot(ux)))
    
        ct = Dd.T.dot(R.dot(ux))
        ev = sly[2:4]-sly[6:8] #velocity error
        #torque to be applied
        #Tau = k1*st*np.sign(ct)
        Tau = k1*st/ct
        #force to be applied
        Phat =R.T.dot(p0-sly[4:6]) #posiciones a ejes cuerpo
        Pdothat =R.T.dot(sly[6:8]) #velocidades a ejes cuerpo
        
        Pr =np.array([Phat[0],Pdothat[0]]) #estados del estimador reducido 
        F = K[0][0:2].dot(Pr)+ K[0][2]*sly[13]
        #F = R.T.dot(k2*ep)[0]-mu.dot(R.T.dot(k3*ev))[0]#
        #F = k2*ct*nD-mu.dot(R.T.dot(ev))[0]# 
        Td[j] = (F + Tau)/2
        Ti[j] = (F - Tau)/2
    pl.figure(5)
    pl.plot(sol.t,Td)
    pl.plot(sol.t,Ti)
    pl.legend(['$T_{S}$','$T_{P}$']) #P port babor izquierda, S starboard
    pl.xlabel('t(s.)')
    pl.ylabel('N.')
    pl.title('$\\theta$ =%2.0f$^o$, r=%2.1fm., v=%2.1fm/s' %(theta*180/np.pi, r, v))
    pl.savefig('./figures/thruster_th%2.0f_r%dm_v%dms.png' %(theta*180/np.pi, r, v))
