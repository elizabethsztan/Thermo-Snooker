#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:31:52 2021

SIMULATION CLASS

@author: liz
"""

class Simulation: 
    
    """
    Simulation takes in 3 arguments
    - num: number of balls
    - R: radius of container
    - r: radius of balls
    """
    
    def __init__ (self, num = 10, R = 10, r = 1): 
    
        self._ball = [Ball(1e20,-R)]#container has big mass so collide affects neglibly
        self._ball = np.array(np.append(self._ball, ball_list(num, r))) #appending balls to ball list
        self._num = num #number of balls, not including container
        #container is first element in self._ball, but not included in self._num count
        
        #for pressure calculation later
        self._impulse = 0 
        self._duration = 0
        
        #for ke/momentum conservation calculation later
        self._intervals = []
        self._total_ke = []
        self._total_momentum = []
        
        #for Maxwell-Boltzmann plot later
        self._total_velocities = []
        
        ball_pos_array = rand_pos_list(num, R, r) #random positions of balls
        ball_vel_array = rand_vel_list(num) #random velocities of balls
        
        for i in range(1,len(self._ball)): #setting pos and vels to balls
            self._ball[i].setpos(ball_pos_array[i-1])
            self._ball[i].setvel(ball_vel_array[i-1])

        #print(self._ball[1]._pos)
        #print(self._ball[1]._vel)

    def next_collision(self):
        
        """
        This functions produces an array of time_to_collision for any ball with each 
        other and the container. Then takes the smallest positive value.
        Moves balls/frame by that value and then collides the two balls/ball and container.
        """
        
        num = len(self._ball)
        time_array = np.zeros((num,num))
        time_array.fill(np.nan)
        
        for i in range(num):
            for j in range(num):
                if i > j:
                    time_array[i][j] = self._ball[i].time_to_collision(self._ball[j])
        #print(time_array)
        tmin = np.nanmin(time_array)
        #print(tmin)
        index = np.where(time_array == tmin) #find which ball collision is soonest
        i = index[0][0] #ball 1/container
        j = index[1][0] #ball 2/container
        
        for ball in self._ball[1:]: #move all balls by tmin
            ball.move(tmin)
        
        #Calculate pressure by finding total impulse of wall on ball, divide by duration
        #to get force, then divide by circumference of container (area) to get pressure
        mom_bef = 0
        mom_aft = 0
        
        if i == 0: #either i or j = 0 for container
            mom_bef = self._ball[j].momentum()
        if j == 0:
            mom_bef = self._ball[i].momentum()

        self._ball[i].collide(self._ball[j]) #collide the two balls
        
        if i == 0:
            mom_aft = self._ball[j].momentum()
        if j == 0:
            mom_aft = self._ball[i].momentum()
        
        #print (mom_bef, mom_aft)
        
        delta_mom = abs(np.dot(mom_aft-mom_bef,mom_aft-mom_bef))
        self._impulse += delta_mom
        
        #find total ke and momentum of the balls after collision
        for i in range(self._num):
            
            total_ke = 0
            total_ke += self._ball[i+1].ke() #total kinetic energy of all balls
            
            total_mom = np.array([0,0])
            total_mom = total_mom + self._ball[i+1].momentum()
            
            self._total_velocities.append(np.sqrt(abs(np.dot(self._ball[i].vel(), self._ball[i].vel()))))
            
        #make list of total ke and momentum of balls - should remain almost constant
        self._total_momentum.append(total_mom)
        self._total_ke.append(total_ke) 
        
        self._intervals.append(self._duration+tmin) #time this collision takes place
        self._duration += tmin #total time
            
        
    
    def run(self, num_frames, animate=True): #this doesn't work without animate=True... why?
        
        if animate:
            f = pl.figure(figsize = [5,5])
            ax = pl.axes(xlim=(self._ball[0].rad,-self._ball[0].rad),ylim=(self._ball[0].rad,-self._ball[0].rad))
            ax.add_patch(self._ball[0].get_patch()) #draw on container
            
            for ball in self._ball[1:]:
                ax.add_patch(ball.get_patch()) #draw on balls
                
        for frame in range(num_frames): #each frame has one collision
            self.next_collision() 
            if animate:
                pl.pause(0.01)
                pl.title(frame)
        if animate:
            pl.show()
    
    def temp(self): #temp depends on the ke 
        
        av_ke = 0
        for i in range(self._num):
            av_ke += self._ball[i+1].ke()
        av_ke /= self._num
        T = (av_ke/1.38)*1e23
        return T
    
    def pressure(self): 
        
        force = self._impulse/self._duration
        p = force/(np.pi*2*-self._ball[0].rad)
        return p
    