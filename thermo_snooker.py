#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 09:31:10 2021

@author: liz

THERMODYNAMICS SNOOKERED

Project Aims:
Write a code to describe a single ball bouncing inside a circular container

Incorporate this into an animation

Extend the code to include multiple balls

Calculate the temperature and pressure

Study the effects of varying the various parameters to deduce the basic laws of thermodynamics
  
"""
#%%
import numpy as np
from numpy.random import randn
import pylab as pl
import math


#Preset exception comments 
vec_exc_2d = "Not right sized vector, only 2-D allowed"
no_col_exc = "The balls don't collide"
wrong_obj_exc = "The object is not of the right class"

#%%

class Ball:
    #make sure to set a container as a ball with neg radius
    def __init__(self, m, rad, pos=[0,0], vel=[0,0]):
        self._pos = np.array(pos, dtype="float64")
        self._vel = np.array(vel, dtype="float64")
        self.rad = rad
        self.m = m
        if self._pos.shape != (2,) or self._vel.shape != (2,):
           raise Exception(vec_exc_2d)
    
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
        return self._pos
    
    def time_to_collision(self, other):
        r = self._pos-other._pos
        v = self._vel-other._vel
        R = self.rad + other.rad
        
        disc = np.dot(r,v)*np.dot(r,v) - (np.dot(v,v))*(np.dot(r,r)-R*R)
        
        time = np.inf #returns nan if balls don't collide
        
        if disc >= 0 and np.dot(v,v) != 0: #take smallest positive real dt
            dt_pos = (-np.dot(r,v)+disc)/(np.dot(v,v))
            dt_neg = (-np.dot(r,v)-disc)/(np.dot(v,v))
        
            if dt_pos<1e-9:
                dt_pos = 0
            if dt_neg<1e-9:
                dt_neg = 0
                
            if self.rad>0 and other.rad>0: #two balls
                if dt_pos>0 and dt_neg>0:
                    time = dt_neg
            
            else: #ball and container
                if dt_pos>0 and dt_neg>0:
                    time = dt_neg
                if dt_pos>0 and dt_neg<1e-9:
                    time = dt_pos
                
            
            
            # if ball on ball
            #   if 2 solutions are phased ie one theyre different, give nan
            #   else use the same collision conditions as before
            #
            # if ball on container
            #   the ball is always phased through the container so ignore the initial conditions
        return time

#Two balls colliding with each other  
    def collide(self, other):
        r1 = self._pos
        v1 = self._vel
        r2 = other._pos
        v2 = other._vel
        m1 = self.m
        m2 = other.m
        
        # print('***')
        # print( r1, v1, r2, v2, m1,m2)
        self._vel = v1 - (2*m2/(m1+m2)) * ((np.dot((v1-v2), (r1-r2)))/(np.dot((r1-r2),(r1-r2))))*(r1-r2)
        other._vel = v2 - (2*m1/(m1+m2)) * ((np.dot((v2-v1), (r2-r1)))/(np.dot((r2-r1),(r2-r1))))*(r2-r1)
        
        #print( r1, v1, r2, v2, m1,m2)
        
    def get_patch(self): #for animation in run method in Simulation class
        if self.rad > 0: #balls have pos radius and are red
            patch = pl.Circle(self.pos(), self.rad, fc="r") 
        else: #container has neg radius and is not filled
            patch = pl.Circle(self.pos(), -self.rad, fc="black", fill = False, ls = "solid")
        return patch
    
    def kinetic_energy(self):
        ke = 0.5*self.m*np.dot(self._vel, self._vel)
        return ke


# Makes list of random velocities
def rand_vel_list(num):
    vels = []
    for i in range (num):
        rand_vel = np.random.uniform(-2,2,2)
        vels.append(rand_vel)
    return vels

def rand_pos_list(num): #make sure they don't overlap
    pos = []
    for i in range (num):
        x = np.random.uniform(-9,9)
        y = np.random.uniform(-np.sqrt(81-x*x),np.sqrt(81-x*x))
        pos.append(np.array([x,y]))
    return pos

def ball_list(num):
    ball_list = []
    for i in range(num):
        ball_list.append(Ball(1,1))
    return ball_list

class Simulation: 
    def __init__ (self, num): #balls is a list of ball objects, vels is list of velocities
        self._ball = [Ball(1e20,-10)]#container
        self._ball = np.array(np.append(self._ball, ball_list(num))) #appending balls to ball list
        self._num = num
        
        ball_pos_array = rand_pos_list(num) #random positions of balls
        print(ball_pos_array)
        ball_vel_array = rand_vel_list(num) #random velocities of balls
        
        for i in range(len(self._ball)-1):
            self._ball[i+1].setpos(ball_pos_array[i])
            self._ball[i+1].setvel(ball_vel_array[i])
        print(self._ball)
    
    def move_balls(self, dt):
        for i in range (self._num):
            self._ball[i+1].move(dt)
    
    def next_collision(self):
        num = self._num+1
        time_array = np.zeros((num,num))
        time_array.fill(np.inf)
        for i in range(num):
            for j in range(num):
                if i > j:
                    time_array[i][j] = self._ball[i].time_to_collision(self._ball[j])
        # print(time_array)
        
        tmin = np.nanmin (time_array)
        index = np.where(time_array == tmin)
        i = index[0][0]
        j = index[1][0]
        
        for ball in range(1,self._num):
            self._ball[ball].move(tmin)   
            
        self._ball[i].collide(self._ball[j])
        
        return time_array
 
    
    # def __str__(self):
    #     return "(Ball position: %s, %s. Ball velocity: %s, %s. Container radius: %s. Container centre: %s.)" % (self._ball._pos[0], self._ball._pos[1], self._ball._vel[0], self._ball._vel[1], self._container.rad, self._container._pos)
    
    #def __repr__(self): fix this later i cba now
       # return "r = (%s, %s), v = (%s, %s)" % (self._pos[0], self._pos[1], self._vel[0], self._vel[1])
       
    def run(self, num_frames, animate=False):
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10,10),ylim=(-10,10))
            for i in range(self._num):
                ax.add_patch(self._ball[i].get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(0.001)
                pl.show()
        
        
        
        # for frame in range(num_frames):
        #     #self.next_collision()
        #     if animate == True:
        #         pl.pause(0.001)
        #         f = pl.figure(figsize = [5,5])
        #         ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
        #         ax.grid()
        #         for i in range (self._num+1):
        #             ax.add_patch(self._ball[i].get_patch())
        #             ax.set_xticks((np.arange(-10,11,1)))
        #             ax.set_yticks((np.arange(-10,11,1)))
        #             pl.show()
    

"""  
- Check whether the random array of balls work
- make a matrix for times between collisions of balls and container, get rid of diagonal
- use matrix to get next_collision to work with multiple balls
- redo the pressure equation, this is not an ideal gas
    
    
"""




