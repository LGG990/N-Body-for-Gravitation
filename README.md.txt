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
