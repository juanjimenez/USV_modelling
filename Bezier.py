#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 16:39:51 2023
Bezier is your friend
This module define functions to obtain Bezier's curves for interpolating 
a set of points.
we impose continuity in the first and second derivatives. in this way, we obtain
smooth curves which can be followed by whatever  autonomous piece of junk
(still in progress :( )
@author: juanjimenez
"""
import numpy as np
import math
import matplotlib.pyplot as pl
def bezier(p,dt):
    """generate a Bezier's curve using the control points
    supplied in p and the parameter dt
    p should be a list of control points (column vectors), and dt the 
    step size of the parameter """
    #calculate first the Binnomal coefficients 
    cbin = []
    u = p.shape[1] #number of control points
    T = np.arange(0,1+dt,dt)
    
    V = np.zeros([u,len(T)])
    for d in range(0,u):
        cbin.append(math.comb(u-1,d)) #comb calculates (u-1)!/(d!(u-1-d)!)
    i = 0
    for t in T:
        for j in range(u):
            V[j,i] = cbin[j]*(1-t)**(u-1-j)*t**j 
            
        i = i+1
        
        cbz = p@V
        #print(t)
        #print(V)
        #print(i)
            
    return cbz,V,T

def bezierpnt(pc,s):
    """returns the value of a Bezier's curve p = B(s) with control points
    defined in pc
    pc shoul be a matrix (2Xn) n = number of control points
    s a parameter value in teh interval [0,1]"""
    cbin=[]
    n = pc.shape[1]
    for d in range(0,n):
        cbin.append(math.comb(n-1,d)) #comb calculates (u-1)!/(d!(u-1-d)!)
    V = np.zeros([n,1])
    for j in range(n):
        V[j] = cbin[j]*(1-s)**(n-1-j)*s**j 
            
        
        
        p = pc@V
        #print(t)
        #print(V)
        #print(i)
            
    return p
    
    
def devbezier(pc,dt=0):
    ''' Returns the control points of a Bezier's curve derivative. If dt>0 it also
    Return the values of the Bezier's curve' for T = [0,1] with step dt'''
    u = pc.shape[1] #number of control points
    n = u-1 #number of control points of its derivative
    pd = np.zeros([2,n]) #empty matrix for the control points of the derivative
    for i in range(n):
        pd[:,i]= (u-1)*(pc[:,i+1]-pc[:,i])
    if dt > 0:
        cbz,V,T = bezier(pd,dt)
        return(pd,cbz,V,T)
    else:
        return pd
        
        
       
    
    


   
# def interpbz(P,dt,dg):
#     """interpolation using bezier curves. P should be a list of point, the USV
#     has to pass through. We employ the same degree for all polynomies (dg)"""
    
#     #first generate a list of control points, based on the data in P
#     #first dimension fit the dimension of the espace P belong to.
#     #second dimension set of control points for each bezier curve
#     L = []
#     for p in range(len(P)):
        
#         L.append()
        
    
def bezierdraw(p,l):
    
    #dibujamos los puntos de paso
    pl.plot(p[0,:],p[1,:],'o-')
    #y la curva de bezier
    pl.plot(l[0,:],l[1,:])      
        
       
        
    
    