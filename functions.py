""" This module contains all the functions for replica exchange script """

# Importing modules
import numpy as np
from numba import njit

# Functions

@njit
def U(x) :
    """ Potential energy """

    a = 4                                                                               # Tune the depth of the well
    U = a * (np.power(x, 2) - 2) * np.power(x, 2)

    return U



@njit
def eval_forces(replicas, dx, v) :
    """ Function to evaluate forces with Langevin thermostat"""
    
    gamma = 0.5                                                                                             # Friction coefficient
    eta = np.random.normal(loc=0.0, scale=1.0, size=replicas.shape[0])                                      # Random number for thermostat noise

    F_cons = - ( U(replicas[:,1] + dx) - U(replicas[:,1] - dx) ) / ( 2 * dx )                               # Conservative forces
    F_friction = - gamma * v                                                                                # Friction
    F_brownian = np.sqrt(2 * gamma * replicas[:,0]) * eta                                                   # Brownian forces
    replicas[:,3] = F_cons + F_friction + F_brownian                                                        # Update forces

    return replicas[:,3]



@njit
def velocity_verlet(replicas, dt, dx, x_max, x_min) :
    """ Velocity-Verlet algorithm """

    v_dt_2 = np.zeros(replicas.shape[0])                                                # Intermediate velocity vector
    
    v_dt_2 = replicas[:,2] + 0.5 * replicas[:,3] * dt                                   # Update velocities pt.1
    replicas[:,1] = replicas[:,1] + v_dt_2 * dt                                         # Update positions
    replicas[:,3] = eval_forces(replicas, dx, v_dt_2)                                   # Update forces
    replicas[:,2] = v_dt_2 + 0.5 * replicas[:,3] * dt                                   # Update velocities pt.2

    return replicas



@njit
def replica_exchange(replicas) :

    # Choose at random one of the possible simulations (except the last one, because if we choose the replica n-1, we try to exchange the n-1 with the n-th replica, which does not exist)
    chosen_replica = np.random.randint(replicas.shape[0] - 1)

    # Compute the probability
    U_i = U(replicas[chosen_replica][1])
    U_f = U(replicas[chosen_replica + 1][1])
    T_i = replicas[chosen_replica][0]
    T_f = replicas[chosen_replica + 1][0]
    
    p = np.exp( - (U_i - U_f) * ( 1/T_i - 1/T_f ) )

    # Try to swap
    if np.random.rand() < p :
        tmp = replicas[chosen_replica + 1][1]
        replicas[chosen_replica + 1][1] = replicas[chosen_replica][1]
        replicas[chosen_replica][1] = tmp

    return replicas
