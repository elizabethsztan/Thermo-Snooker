#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 08:56:58 2021

BALL CLASS

@author: liz
"""

#%%
import numpy as np
#from numpy.random import randn
import matplotlib.pyplot as plt
import pylab as pl
import math
from scipy import optimize

#Preset exception comments 
vec_exc_2d = "Not right sized vector, only 2-D allowed"
no_col_exc = "The balls don't collide"
wrong_obj_exc = "The object is not of the right class"

#%%

class Ball:
    
    def __init__(self, m, rad, pos=[0,0], vel=[0,0]): #make sure to set a container as a ball with neg radius
        
        self._pos = np.array(pos, dtype="float64")
        self._vel = np.array(vel, dtype="float64")
        self.rad = rad
        self.m = m
        self.initialpatch = None
        if self._pos.shape != (2,) or self._vel.shape != (2,):
           raise Exception(vec_exc_2d) #to make sure pos, vel are of the right dimensions
    
    def __str__(self):
        return "(Position: %s, %s. Velocity: %s, %s)" % (self._pos[0], self._pos[1], self._vel[0], self._vel[1])
    
    def __repr__(self):
        return "r = (%s, %s), v = (%s, %s)" % (self._pos[0], self._pos[1], self._vel[0], self._vel[1])
    
    def pos(self):
        return self._pos
    
    def vel(self):
        return self._vel
    
    def setpos(self, newpos):
        self._pos = (newpos)
        if self._pos.shape != (2,):
            raise Exception(vec_exc_2d)
    
    def setvel(self, newvel):
        self._vel = (newvel)
        if self._vel.shape != (2,):
            raise Exception(vec_exc_2d)
    
    def move(self, dt):
        self._pos = self._pos + self._vel*dt
        try:
            self._patch.center = self._pos
        except:
            print('nop')
        return self._pos
    
    def time_to_collision(self, other): #BIG BIG BIG 
        
        r = self._pos - other._pos 
        v = self._vel - other._vel
        R = self.rad + other.rad
        
        #dot products for calc
        dot_rr = np.dot(r,r)
        dot_rv = np.dot(r,v)
        dot_vv = np.dot(v,v)
        
        disc = dot_rv*dot_rv-dot_vv*(dot_rr-R*R)
        
        if disc >= 0 and dot_vv>0: #get the real values 
            dt_pos = ((-dot_rv) + np.sqrt(disc))/dot_vv
            dt_neg = ((-dot_rv) - np.sqrt(disc))/dot_vv
        
            if dt_pos>0 and dt_neg >0: #only pos time
            
                if self.rad<0 or other.rad<0: #ball and container collision
                    time=dt_pos 
                else: #ball and ball collision
                    time=dt_neg
                    
            else:
                if dt_pos<0 and dt_neg<0: #no future collision
                    time=np.NaN 
                else:
                    if self.rad<0 or other.rad<0: #ball and container collision
                        time=dt_pos
                        
                    else: #when balls overlap, returns nan
                        time=np.NaN
            return time
        
    #collision of two balls
    def collide(self, other):
        
        r1 = self._pos
        v1 = self._vel
        r2 = other._pos
        v2 = other._vel
        m1 = self.m
        m2 = other.m
        
        self._vel = v1 - (2*m2/(m1+m2)) * ((np.dot((v1-v2), (r1-r2)))/(np.dot((r1-r2),(r1-r2))))*(r1-r2)
        other._vel = v2 - (2*m1/(m1+m2)) * ((np.dot((v2-v1), (r2-r1)))/(np.dot((r2-r1),(r2-r1))))*(r2-r1)
    
    def ke (self):
        ke = 0.5*self.m*np.dot(self._vel,self._vel)
        return ke
    
    def momentum (self):
        p = self.m*self._vel
        return p
    
    def get_patch(self): #for animation in run method in Simulation class
    
        if self.rad > 0: #balls have pos radius and are red
            self._patch = pl.Circle(self._pos, self.rad, fc="r") 
        else: #container has neg radius and is not filled
            self._patch = pl.Circle(self._pos, -self.rad, fc="black", fill = False)
        return self._patch
    
    
    
    
    
    