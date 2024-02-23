#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 15:48:13 2023
This module define several functions to use gvfparamatric wit Bezier's curves
@author: juanjimenez
"""
import Bezier as Bz
import numpy as np
import math
import matplotlib.pyplot as pl

def bzxi(x,pt,s,k):
    ''' Calculates de augmented vector field (3D) from the derivative of 
    a Bezier'curve that defines the desired trayectory
    x = point to calculate the value of the field
    pt Bezier's curve control points
    s Bezier's curve parameter s in (0,1)
    k gain for the normal component to the curve of the vector field
    
    ''' 
    if len(x.shape)==3:
        
        #print(x)
        pd = Bz.devbezier(pt) #control points of Bezier's curve derivative
        p  = Bz.bezierpnt(pt,s) #value of f(s)
        dp = Bz.bezierpnt(pd,s) # value of f'(s) 
        phi = np.transpose(x,(1,0,2))-p #value of phi (surfaces to define de bezier's curve)
        #print(phi)
        #grd gradient of phi1 and gradient of phi2, both puth togueder in a 3x2
        #matri grd = [grad(phi1),grad(phi2)]
        kph = k*phi
        #print(kph)
        grd = np.vstack(([1,0],[0,1],[-dp[0,0],-dp[1,0]])) 
        xi = (np.vstack((dp,1)).T-grd.dot(kph).T).T
        #print(xi)
    
        return xi
    else:
        pd = Bz.devbezier(pt) #control points of Bezier's curve derivative
        p  = Bz.bezierpnt(pt,s) #value of f(s)
        dp = Bz.bezierpnt(pd,s) # value of f'(s) 
        phi = x-p #value of phi (surfaces to define de bezier's curve)
        #grd gradient of phi1 and gradient of phi2, both puth togueder in a 3x2
        #matri grd = [grad(phi1),grad(phi2)]
        kph = k*phi
        grd = np.vstack(([1,0],[0,1],[dp[0,0],dp[1,0]])) 
        xi = np.vstack((dp,1))-grd.dot(kph)
    
        return xi
    
    