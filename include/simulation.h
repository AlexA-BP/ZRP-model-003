#ifndef SIMULATION_H_785875888084
#define SIMULATION_H_785875888084

#include <stddef.h>
#include "pcg_variants.h"

void sim(pcg32_random_t* rngptr, int* particles, int particles_size, int* lat, 
        int lat_size, double dt, double b, double gamma); 

void init_particles(pcg32_random_t* rng, int *arr, int L, int N, int S);

void init_lattice(int *particles, int particles_size, int *lat, int lat_size, int species); 

void sim_test(pcg32_random_t* rngptr, int* particles, int particles_size, int* lat, 
        int lat_size, double dt, double b);
    
#endif