#!/usr/bin/env python
# coding: utf-8

# In[5]:


###### Program to simulate the gravitational motion of bodies, given their initial conditions of: mass,position, and velocity.
#The data is read in from an external file using np.loadtxt from the numpy library.
#Values are read as a string and converted to the correct data types as required.
#The motion of the bodies is calculated using the Velocity Verlet method. Given a maximum simulation time and timestep,
#the program produces a plot of position X vs Y of each body at every timestep, a plot of the kinetic, potential and
#total energy vs time, reproduces Keplers third law: a^3/T^2 = const. for each body, prints the angular momentum at the
#final time step, and the error in the total energy.

#Three files are included for simulation: fig8.txt this is the stable 3 body ; earth.txt earth sun orbit ; solar system.txt

#Recommended simulation times:

#fig8.txt   TIME = >3 seconds, STEP = <0.01 seconds 
#earth.txt  TIME = >31540000 seconds (a calendar year), STEP = 86400 seconds (a day)
#solar system.txt Time = >930000000 seconds (Saturn orbital period), STEP = 86400 seconds

#imports
import numpy as np
import math
import matplotlib.pyplot as plt

#which file to read, simulation time, and timestep from user input
fileName = input("Enter file name, including extension\n")
timeMax = int(input("Enter the number of seconds to simulate to as an integer\n"))
step = float(input("Enter the step size\n"))

#Read data file
data = np.loadtxt(fileName,skiprows=1,delimiter = '\t',dtype='str')

#user input which body to simulate energies
print("There are ",len(data)," bodies loaded\n")
simulateE = int(input("\nWhich # body would you like to calculate the energy and angular momentum for?\n"))-1


#Define constants
time = np.linspace(0,timeMax,int((timeMax/step)))
numB = len(data)
sizeT = len(time)

#G is different for fig 8 stable orbit simulation
if fileName == 'fig8.txt':
    G = 1
else:
    G = 6.67*10**(-11)

#Define variables
mass = np.zeros((numB))
name = []
position = np.zeros((sizeT,numB,3)) #indexing: time, particle, xyz element
velocity = np.zeros((sizeT,numB,3))
velocity_half = np.zeros((sizeT,numB,3))
acceleration = np.zeros((sizeT,numB,3))
r = np.zeros((sizeT,numB,numB,3)) # indexing: time, ith particle, distance between ith and jth particle, xyz

#Array r is created such that for every time step the distances from each partcle are caluclated from every perspective
#The gravitational force acting on body A is calculated using the distances AB and AC but not BC. Having the array
#in this setup makes the association of each distance for each gravitational acceleration caluclation much easier

#Populate variables from file
for i in range(numB):
    name.append(data[i][0])
    mass[i] = data[i][1].astype('float')
    
    for element in range(3):
        position[0][i][element] = data[i][element+2].astype('float')
        velocity[0][i][element] = data[i][element+5].astype('float')

def sumSqCubed(vector): # given a vector, this function squares the elements, sums them, and returns
                        # the 3/2 power of the result
    temp = np.zeros((1))
    for i in range(3): #loop elements
        temp[0] += float(vector[i]**2)
    return(temp[0]**(3/2))
    
def sumSqrt(vector): #given a vector, this function squares the elements, sums them, and returns
                     # the 1/2 power of the result
    temp = np.zeros((1))
    for i in range(3):#loop elements
        temp[0] += float(vector[i]**2)
    return(temp[0]**.5)     


#rValue function, given the current position of every particle and the number of particles, calculates the distances
#between particles. All the distances are only calculated once otherwise distance AB will be caluclated twice when 
#distance BA is calculated. Distances from a particle to itself are also neglected.
# Makes use of the nested summation 
#                                        Sum_j_from_1_to_n(  Sum_i_from_j+1_to_n(  P_j - Q_i ) )

#Note: the variables P and Q represent the position vectors, where n is the maximum number of particles
def rValue(position,numB):  
    r = np.zeros((numB,numB,3))
    for i in range(numB): # loop ith particle
        for j in range(numB): # loop jth particle
            if i < j: #avoid calculating distance from itself
                for x in range(3): #loop elements
                    r[i][j][x] = position[i][x] - position[j][x]
                    r[j][i][x] = -r[i][j][x]
    return r

def Acceleration(mass,timestep,r,G):#given mass, the current timestep, disctances, and G
                                  # calculates gravitational acceleration of particle i due to all other particles
                                  # acceleration is returned as a vector after summing x,y,z contributions 
                                  # from j particles
                                  # a = -Gm r_vector/r^3
                    
    numB = len(mass)
    accel = np.zeros((numB,3))
    
    for i in range(numB): # ith particle
        for n in range(numB): # nth particle
            for e in range(3):
                if n != i: # make sure not to use its own mass in calculation
                    accel[i][e] -= G*mass[n]*r[i][n][e]/sumSqCubed(r[i][n])
    return(accel) #returns acceleration of all particles, at the given time, as a vector

###Energy and angular momentum calculation
def Energy(mass,position,velocity,numB,t,body,G):
    totalMass = np.zeros((1))
    Sum = np.zeros((3))
    Pe = np.zeros((numB))
    Ke = np.zeros((numB))
    totalE = np.zeros((numB))
    angMom = np.zeros((numB,3))
    com = np.zeros((3))
    
    #calculate center of mass
    for i in range(numB): #loop particles
        totalMass += mass[i]
        for x in range(3): # loop elements
            Sum[x] += position[t][i][x] * mass[i]     
    com = (1/totalMass)*Sum
    
    #calculate potential and kinetic enrgy
    for i in range(numB): #loop ith particles
        Ke[i] = .5*mass[i]*sumSqrt(velocity[t][i])**2
        for body in range(numB): #loop all particles and avoid ith particle
            if i != body:
                r = position[t][i] - position[t][body]
                Pe[i] -= G*mass[i]*mass[body]/sumSqrt(r)

    for i in range(numB): #loop particles
        totalE[i] = Pe[i] + Ke[i]
        angMom[i] = abs(position[t][i]- com) * mass[i]*velocity[t][i]
    
    return totalE[i],Pe[i],Ke[i],angMom[i]

#Calculates semi major axis and time period using Keplers law. Returns a^3/T^2, expected to be constant for all
#bodies orbiting a common center of mass
def thirdLaw(position,body,mass,name,G):
    posBody = np.zeros((len(position)))
    for t in range(len(position)): #grabs all x positions of the chosen body
        posBody[t] = sumSqrt(position[t][body] - position[t][0]) #finds the max distance between the body and the sun
    a = abs(max(posBody)) # finds the maximum value of position as the semi major axis
    T = math.sqrt(((a**3)*4*(math.pi)**2/(G*mass[0])))
    if a and T != 0:
        print("\n\nRatio of cube of semi-major axis with square of time period for ",name[body]," is ", a**3/T**2)
    return

### Velocity Verlet calculation

#calculate initial step to kick the calculation loop
r[0] = rValue(position[0],numB)
acceleration[0] = Acceleration(mass,0,r[0],G)
velocity_half[0] = velocity[0] + 0.5*step*acceleration[0]

for t in range(1,len(time)): # loop time step
                             # for every step calcualate new values of velocity half step, position then acceleration
                             # using the acceleration function, then the new velocity
    velocity_half[t] = velocity[t-1] + .5*step*acceleration[t-1]
    position[t] = position[t-1] + step*velocity_half[t]
    r[t] = rValue(position[t],numB)
    acceleration[t] = Acceleration(mass,t,r[t],G)
    velocity[t] = velocity_half[t] + .5*step*acceleration[t]

###Plotting
plt.figure(figsize=(14,14))

bodies = np.zeros((numB,2,sizeT)) # array to essentially change the shape of position in order to slice over all
                                  # time steps when plotting
for t in range(len(time)):
    for n in range(numB):
        for element in range(2):
            bodies[n][element][t]=(position[t][n][element])
            
for n in range(numB): # plot all bodies. Scatter plot the sun for orbital simulations as it doesn't show up
                      # well in plt.plot because it barely moves in comparison
    if n == 0:
        if fileName == 'fig8.txt':
            plt.plot(bodies[0][0],bodies[0][1],label=name[n])
        else:
            plt.scatter(bodies[0][0],bodies[0][1],label=name[n])
    else:
        plt.plot(bodies[n][0],bodies[n][1],label=name[n])

plt.xlabel('X [m]',fontsize = 22)
plt.ylabel('Y [m]',fontsize = 22)
plt.legend(fontsize=18)
plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20)
plt.show()

##Energy vs Time data popoulated from the Energy function
totE = np.zeros(sizeT)
Ke = np.zeros(sizeT)
Pe = np.zeros(sizeT)
angMom = np.zeros((sizeT,3))
for t in range(len(time)):
    totE[t],Pe[t],Ke[t],angMom[t] = Energy(mass,position,velocity,numB,t,simulateE,G)   
    
plt.figure(figsize=(14,14))
plt.plot(time,totE,label='Total Energy')
plt.plot(time,Pe,label='Potential Energy')
plt.plot(time,Ke,label='Kineic Energy')
plt.legend(fontsize=16)
plt.title('Energy Vs Time',fontsize=25)
plt.ylabel('Energy [ J]',fontsize=22)
plt.xlabel('Time [S]',fontsize=22)
plt.xticks(fontsize=20,rotation=45)
plt.yticks(fontsize=20)
plt.show()

#Keplers third law and angular momentum
if fileName != 'fig8.txt':
    for i in range(numB):
        thirdLaw(position,i,mass,name,G)


#print angular momentum and energy calculation error
#error is calculated using the total energy: initial-final/final *100

final = int(timeMax/step)-1
print("\n\nAngular Momentum of ",name[simulateE]," after ",timeMax," seconds =\n\n ",angMom[final],"[Kg m^2 s^-1]")
print("\n\n\nError in simulation with timestep of ",step," is \n", abs(abs(totE[1]-totE[final])/totE[final])*100,"%")



#Future additions:#########################################
#Create 3-D plot
#Animate the 3-D plot over the whole simulation time
#Fix the potential energy calculation, currently calculates twice I believe

