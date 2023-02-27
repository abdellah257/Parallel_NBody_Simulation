# Parallel Simulation of N-Body System

## Description:

This is an implementation of an MPI version of Barnes Hut Algorithm for N-Body systems. Here is a description of the project folder :

| File                  |                                               Description | Inputs                                    |
| :-------------------- | --------------------------------------------------------: | :-----------------------------------------|
| **lib.h**             |    Contains structures definition used for the simulation |                                           |
| **bh_sim_seq.cpp**    |                           Sequential Barnes-Hut Algorithm | arg[1]: N_particles, arg[2]: N_iterations |
| **bh_sim_par.cpp**    |                             Parallel Barnes-Hut Algorithm | arg[1]: N_particles, arg[2]: N_iterations |

The files speed_seq.txt and speed_par.txt contain the execution time for diferent numbre of particles.

## Compilation

You can compile and execute the MPI program using :
```bash
$ mpic++ bh_sim_par.cpp -o bh_sim_par
$ mpirun -n <number_of_nodes> bh_sim_par <Numbre of particles> <Numbre of iterations>
```

You can compile and execute the sequential program using :
```bash
$ g++ bh_sim_seq.cpp -o bh_sim_seq
$ ./bh_sim_seq <Numbre of particles> <Numbre of iterations>