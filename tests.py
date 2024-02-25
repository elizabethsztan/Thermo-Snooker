#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 15:02:38 2021


TESTS


@author: liz
"""


#Quickdraw functions to test location of ball/quickly animate

#Draw one ball and container
def quickdraw1(ball1, container):
        f = pl.figure(figsize = [5,5])
        patch1 = pl.Circle(ball1.pos(), ball1.rad, fc="r")
        patchc = pl.Circle(container.pos(), -container.rad, fc="black", fill = False, ls = "solid")

        ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))

        ax.add_patch(patch1)
        ax.add_patch(patchc)
        ax.grid()

        ax.set_xticks((np.arange(-10,11,1)))
        ax.set_yticks((np.arange(-10,11,1)))
        
#Draw on position of two balls 
def quickdraw2(container, ball1,ball2):
        f = pl.figure(figsize = [5,5])
        patch1 = pl.Circle(ball1.pos(), ball1.rad, fc="r")
        patch2 = pl.Circle(ball2.pos(), ball2.rad, fc="b")
        patchc = pl.Circle(container.pos(), -container.rad, fc="black", fill = False, ls = "solid")

        ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))

        ax.add_patch(patch1)
        ax.add_patch(patch2)
        ax.add_patch(patchc)
        ax.grid()

        ax.set_xticks((np.arange(-10,11,1)))
        ax.set_yticks((np.arange(-10,11,1)))


#Draw on positions of balls moving for set time w interval dt (integer dt please)
def animate(ball1, ball2, container, dt, time):
    for i in range (dt,time):
        ball1.move(dt)
        ball2.move(dt)
        quickdraw(ball1,ball2,container)



#%%


sim3 = Simulation(Ball(1e20,-10), 3)
sim3.run(1,True)

sim2._ball[0].setpos([-5,0])
sim2._ball[0].setvel([1,0])

sim2._ball[1].setpos([3,8])
sim2._ball[1].setvel([-2,-5])
#%%

sim = Simulation (1)

sim.run(1, True)
print(sim.next_collision())

#%%
sim.move_balls(0.78873568328452)
sim.run(1, True)




#%%
num = 10
time_array = np.zeros((num,num))
time_array.fill(np.NaN)
for i in range(num):
    for j in range(num):
        if i > j:
            time_array[i][j] = ball[i].time_to_collision(ball[j])
print(time_array)


#%%

#Testing the next_collision() method

ball = Ball(1,1)
container = Ball(1e20,-10)

sim1 = Simulation (container, ball)

sim1._ball.setpos([-5,0])
sim1._ball.setvel([1,0])
print(sim1.pressure())
print (sim1._ball.kinetic_energy())

sim1.run(1,animate = True)

print (sim1._ball.kinetic_energy())

print(sim1.pressure())


"""
quickdraw1(ball, container)

sim1.next_collision()

quickdraw1(ball, container)

sim1.next_collision()

quickdraw1(ball, container)
"""

        
 

#%%



ball1 = Ball(1,1)
ball2 = Ball(1,1)
container = Ball(1e10,-10) #container is ball w neg radius

ball1.setpos([0,4])
ball1.setvel([0,-1])
ball2.setpos([0,-3])
ball2.setvel([0,2])

quickdraw(ball1,ball2, container)

time = ball1.time_to_collision(ball2)

ball1.move(time)
ball2.move(time)

quickdraw(ball1,ball2, container)



ball1.collide(ball2)

animate(ball1,ball2,container, 3)

print (ball1.vel, ball2.vel)



#%%

