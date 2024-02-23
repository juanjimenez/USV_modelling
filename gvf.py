# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 20:30:52 2022
guiado mediante el algoritmo de gvf
@author: abierto
"""
import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as pl
def elipse(a,b,alpha,p,p0):
    pe = p-p0
    R = np.array([[np.cos(alpha),-np.sin(alpha)],[np.sin(alpha),np.cos(alpha)]])
    peb = R.dot(pe)
    e = peb[0][0]**2/a**2+peb[1][0]**2/b**2-1
    n = np.array([[2*peb[0][0]/a**2],[2*peb[1][0]/b**2]])
    H = np.array([[2/a,0],[0,2/b]])
    return e,n,H
    
def circulo(p,p0,r):
    pe = p-p0
    
    e = pe[0][0]**2+pe[1][0]**2-r**2
    n = np.array([[2*pe[0][0]],[2*pe[1][0]]])
    H = np.array([[2,0],[0,2]])
    return e,n,H
def gvf_control_boat(p,dot_p,e,n,H,ke,kd,s,m):
    """
    Parameters
    ----------
    p : 2D vector
        USV position.
    dotp : 2D vector
        USV velocity.
    e : TYPE
        Trayectory to be tracked as a implicit function f(x,y) = e.
        f(x,y) = 0 is the desired trayectory
    n : 2D vector
        gradient of e in p
    H : 2x2 Matrix (defined as a numpy array anyway)
        Hessian of e in p.
    ke : constant
        gain of the gvf in the gradient direction of e.
    kd : constant
       gain to control de convergence rate to the gvf.
    s : constant 
        USV speed IN HEADING DIRECTION, this variable depends only on the USV
        Dynamics and it is supposed that can be always fit by the USV inner
        control loop...
    m: 2D vector
        USV attitude  

    Returns
    -------
    phi_dot : chage rate of the USV heading (yaw) angle It is an input to the 
    usv control sistem...
    """

    #Matriz de rotación noventa grados
    E =np.array([[0,1],[-1,0]])
    #Vector de heading del USV
    
    
    dot_pd = E.dot(n)-ke*e*n
    
    ddot_pd = (E - ke*e*np.eye(2)).dot(H).dot(dot_p) - ke * n.T.dot(dot_p) * n
    ddot_pdhat = -E.dot(dot_pd.dot(dot_pd.T)).dot(E).dot(ddot_pd)/npl.norm(dot_pd)**3
    dot_Xd = ddot_pdhat.T.dot(E).dot(dot_pd)/npl.norm(dot_pd)
    aux = kd*dot_p.T.dot(E).dot(dot_pd)/(npl.norm(dot_p)*npl.norm(dot_pd))
    
    phi_dot = dot_Xd+kd*dot_p.T.dot(E).dot(dot_pd)/(npl.norm(dot_p)*npl.norm(dot_pd))
    
    #Esta es la parte que no funciona y por eso la tengo comentada
    #phi_dot = phi_dot*npl.norm(dot_p)**2/s/dot_p.T.dot(m)
    
    
    
    return phi_dot, dot_Xd, aux
   
    
    
    
    
    
    