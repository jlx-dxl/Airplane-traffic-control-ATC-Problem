# Overview
Airplane traffic control (ATC) is responsible for helping aircraft navigation by identifying routes such that aircraft can get to their destinations as quickly as possible, but at the same time avoiding or preventing collisions with other aircraft.

In this project, you will implement an ATC strategy in which you are given a configuration file that specifies the following for a set of _p_ planes:
* the origin coordinates of the plane's flight
* the destination coordinates of the flight
* the departure time of the flight

The origin and destination are given as Cartesian x/y-coordinates, with (0, 0) as the upper left corner and (100,100) as the bottom right corner; planes may not leave this legal airspace. 
The departure time is given as the number of "steps" after the start of the simulation; a flight may not take off before its departure time, but it may leave afterward.

Two planes that are in the air may not come with 5 units of each other, and planes may not move on the z-axis (just go with it, okay?). 
If two planes in the air come within 5 units of each other, the simulation will terminate.

In each step of the simulation, your ATC strategy will indicate the bearing (direction) of each of the _p_ planes. 
At the start of the simulation, all planes are assumed to be on the ground, and have a bearing of -1. 
When a plane takes off, its bearing changes to its initial direction, with 0 being due north, 90 being due east, 180 being due south, and 270 being due west; a bearing of 360 is allowed (also due north), but after takeoff, any bearing less than 0 or greater than 360 will be considered illegal and the simulation will terminate.

Planes that are in the air move at a constant velocity of one unit per time step. 
Between each step, a plane's bearing may change by a maximum of ±10 degrees: a plane that is going due east (bearing 90) can change to a bearing between 80 and 100, but cannot simply make a sharp right turn and head directly to bearing 180 in one step.

Once a plane reaches its destination, it lands and has its bearing changed to -2 to indicate that it is on the ground. 
The simulation ends when each plane has reached its destination.

**The primary goal is to get all planes to their destinations in the fewest number of steps**, without any planes taking off before their departure time or flying too close to each other.

The amount of time the planes spend in the air is considered to be the amount of "power" used: you should be attempting to minimize power in addition to minimizing the time it takes to get all planes to their destinations. 
Additionally, if a plane stays on the ground after its departure time, this is considered a "delay"; you should attempt to minimize this value, too, but focus on getting all planes to their destinations as quickly as possible.

# Setting Up
One member of your team should [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) this repo so that your team can submit your solution.
Other members of the team should work with that same fork.

Then, [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) your forked repo (or just download the code) to your local computer so that you can work with the code and run the simulator.

Once you have downloaded the code, open IntelliJ IDEA, select "Open" to open an existing project, and select the folder containing the files in this repo.

IntelliJ should compile and build the code for you.

Run "airplane.sim.GameEngine.main()" and you should see the simulator UI launch. 
You can configure:
* the Java class to use for your ATC strategy
* the configuration file indicating each plane’s initial location, destination, and departure time (described below)
* the delay of refreshing the UI

To run a simulation:
* Click "Begin New Game" to start the simulation.
* Click "Step" to move ahead one unit of time.
* Click "Play" to move continuously
* Click "Pause" to temporarily pause the game.
* Click "Resign" to quit.

When the simulation ends, you will see the total simulation time, the amount of power used by all airplanes (time they were in the air), and total delay.

# Implementing Your Solution
Create a class called airplane.gX.GroupXPlayer where _X_ is your group number. 
This class must:
* extend airplane.sim.Player
* implement the getName, startNewGame, and updatePlanes methods

See airplane.g0.SerializedPlayer for an example.
  
The **startNewGame** method is called at the beginning of the simulation, and provides a List of Plane objects, which include the current location (which, at this point, is its origin), destination, current bearing (which will be -1 since it’s on the ground) and departure time of each flight.

Then, at each step in the simulation, the **updatePlanes** method is called. 
The arguments are the List of Plane objects (with their current location and bearing), the round number, and an array of bearings; your method should update and then return that array.

To determine the bearing needed to get from point A to point B, you can call the **calculateBearing** method in the Player superclass.

Your player may also run simulations within the simulation, e.g. to determine whether planes will collide or how close they may get to each other before actually committing to those moves. 
At any point in your player’s execution, it may call the **startSimulation** method. 
This will repeatedly call your player’s **simulateUpdate** method and then update the simulated planes accordingly. 
You may terminate the simulation at any point by calling **stopSimulation** 
Otherwise, the simulation will run until all planes have reached their destinations, at which point startSimulation will finish. 
The return value of startSimulation is a SimulationResult object that indicates the time at which the simulation stopped, the reason for stopping, and the List of Planes at the point when the simulation finished.

To add your Player to the application, add the name of your class to **airplane.xml** in the "airplane .classes" entry (if you have more than one, the class names should be separated by whitespace).

Note: To do logging/debugging, do _not_ use System.out.println or System.err.println. 
Rather, create a Logger instance (see SerializedPlayer for an example) and then call its trace, debug, info, warn, or error method and pass the String to appear in the console.

# Configuration Files
The "flights" directory of the IntelliJ project contains configuration files for different situations that your ATC strategy should be able to address.

Each line of the configuration file contains three semicolon-separated fields, each describing an individual flight:
* the x- and y-coordinates of the initial location (origin)
* the x- and y-coordinates of the destination
* the departure time; the flight is not allowed to leave the origin before this time

We will begin by using these configurations for our initial evaluation of strategies (probably in this order):
* Simple.txt: a single flight
* Parallel.txt: two flights that travel parallel to each other, thus have no chance of collision
* Double.txt: two flights that travel directly toward each other
* Cross.txt: two flights whose paths cross
* Three.txt: three flights all heading to the same destination
* Four.txt: four flights all heading to the same destination

More configurations will be added as we progress through this project, and you are welcome and encouraged to create your own if you think of interesting configurations or situations that your strategy is particularly adept at handling.
