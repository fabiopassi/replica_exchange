# replica_exchange

This repo contains the necessary files to perform the simulation of a 1D replica exchange on a single particle.

## Physics background

This simulation involves a particle moving in a 1D double-well potential.

The chosen ensemble is the NVT and the temperature is fixed with a Langevin thermostat.

The basic idea behind the replica exchange method is the following: usually, we are interested in simulating a certain system at a given temperature $T$, e.g. a protein at 310 °C; this system may have several stable states -separated by energy barriers much higher than the typical thermal energy $ k_B T $- and waiting for a normal molecular dynamics simulation to explore all the states requires definitely too many time, since the probability of overcoming a barrier with height $\Delta$ is proportional to $\exp{-\Delta / k_B T}$.  
> Just to clarify, an example of "several stable states" is given by the many different conformations that a protein.

Nonetheless, one could imagine to simulate the very same system at a much higher temperature $\tilde{T}$: in this way, the barriers could be easily overcome, but the price to pay is that we do not sample the canonical ensemble at temperature $T$, which is our aim.

The solution consists in using some high-temperature copies of the system to easily overcome the barriers, hence generating configurations belonging to different stable states; then, these configurations must be "transferred" to the system at temperature $T$, so that we can sample all the energy minima with the canonical distribution at the correct temperature.

The practical scheme to do so is the following:

1. $M$ copies of the system are created, each one with a different temperature $T_i$ and with $T_M > T_{M-1} > ... > T_1$. These copies are known as "replicas".

2. Each copy is simulated via molecular dynamics or Monte Carlo.

3. Every $N_{exchange}$ steps, we select two neighbouring copies (i.e. with temperatures $T_K$ and $T_{K+1}$) and we try to exchange their coordinates using a probabilistic term, which depends on the value of the potential in the two replicas and on their temperatures. More details can be found in chapter 7.5 of the book "Statistical Mechanics: Theory and Molecular Simulation" by Mark Tuckerman.

> To be more precise, the replica exchange method where the temperature is used as control parameter is known as "parallel tempering".

## Technical details: python version

The scripts are written in python. The only additional packages required to run the simulation are `numpy`, `numba` and `matplotlib`; if you have conda installed, the command:

```bash
conda create -n replica_exchange numba matplotlib
```

should create an environment with all the necessary packages.

After this, you can start the simulation with the commands:

```bash
conda activate replica_exchange
python main.py
```

## Contact

If you have any doubt or if you spot mistakes/bugs, feel free to contact [fabio.passi24@gmail.com](fabio.passi24@gmail.com)
