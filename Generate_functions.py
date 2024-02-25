#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:23:10 2021

GENERATE BALL LIST, POS AND VEL LIST

@author: liz
"""

# Makes list of random velocities
def rand_vel_list(num):
    vels = []
    for i in range (num):
        rand_vel = np.random.uniform(-3,3,2) #3m/s is the max vel in x or y
        vels.append(rand_vel)
    return vels

#Makes list of random positions inside container
def rand_pos_list(num, R, r): #R is radius of container, r is radius of ball
    r_new = R-r
    pos = []
    for i in range (num):
        x = np.random.uniform(-r_new, r_new)
        y = np.random.uniform(-np.sqrt(r_new*r_new-x*x),np.sqrt(r_new*r_new-x*x))
        pos.append(np.array([x,y]))
    return pos
#Note: they have the potential to overlap but this is phased out in time_to_collision

#Generate a list of ball objects
def ball_list(num, r): #r is radius of ball
    ball_list = []
    for i in range(num):
        ball_list.append(Ball(1,r))
    return ball_list