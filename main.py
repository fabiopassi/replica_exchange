##########################################################################################
#####                        REPLICA EXCHANGE 1D (NVT)                                #### 
##########################################################################################

# Importing modules
from functions import *
import numpy as np
import matplotlib.pyplot as plt

# Variables
x_max = 1.6                                                             # Right border
x_min = -x_max                                                          # Left border
dx = 0.0001                                                             # Spacing to evaluate the force
dt = 0.002                                                              # Timestep for Verlet algorithm
t_exchange = 50                                                         # Number of steps every which one performs replica exchange
num_steps = 50000                                                       # Number of steps               
num_temperatures = 20                                                   # Number of temperature values (i.e. number of simulations)
T = 50                                                                  # Temperature at which the system is really simulated
T_max = 1000                                                            # Maximum possible value of temperature
dT = (T_max - T) / (num_temperatures - 1)

replicas = np.zeros((num_temperatures, 4))                              # All the replicas: each row is a replica at different temperature

# System initialization
replicas[:,0] = np.arange(num_temperatures) * dT + T                    # Col 0 : temperatures
replicas[:,1] = np.ones(num_temperatures)                               # Col 1 : positions
replicas[:,2] = np.random.rand(num_temperatures)                        # Col 2 : velocities
replicas[:,3] = np.zeros(num_temperatures)                              # Col 3 : accelerations
print("Temperatures :\n", replicas[:,0])

x_repl_0 = []                                                           # Trajectory of the replica at temperature T
x_repl_1 = []                                                           # Trajectory of the replica at temperature T + dT

# Evaluate initial forces
replicas[:,3] = eval_forces(replicas, dx, replicas[:,2])


# Time evolution
for i in range(num_steps):

    # Time evolution with velocity-verlet algorithm
    replicas = velocity_verlet(replicas, dt, dx, x_max, x_min)

    # Try to exchange neighbouring replicas
    if i % t_exchange == 0 :
        replicas = replica_exchange(replicas)

    # Store position of the simulation at temperature T
    if i % 2 == 0:
        x_repl_0.append(replicas[0][1])
        x_repl_1.append(replicas[1][1])


# Plot
dx_plot = 0.01
x_axis = np.arange(x_min, x_max, dx_plot)
U_plot = U(x_axis)

# to run GUI event loop
plt.ion()
 
# here we are creating sub plots
figure, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 13))
ax1.plot(x_axis, U_plot, color = "blue")
ax2.plot(x_axis, U_plot, color = "blue")
curr_1 = ax1.scatter(x_repl_1[0], U(x_repl_1[0]), color="red")
curr_0 = ax2.scatter(x_repl_0[0], U(x_repl_0[0]), color="orange")

ax1.set_title("Replica T+dT", fontsize=18)
ax2.set_title("Replica T", fontsize=18)
 
# setting x-axis label and y-axis label
figure.supxlabel("x", fontsize=18)
figure.supylabel("U", fontsize=18)
 
# Loop
step_plot = 10
for i in range(0, len(x_repl_0), step_plot):

    # updating data values
    curr_1.remove()
    curr_1 = ax1.scatter(x_repl_1[i], U(x_repl_1[i]), color="red")
    curr_0.remove()
    curr_0 = ax2.scatter(x_repl_0[i], U(x_repl_0[i]), color="orange")
 
    figure.canvas.flush_events()

