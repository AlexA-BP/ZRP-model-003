#include <math.h>

#include "simulation.h"
#include "pcg_variants.h"
#include "entropy.h"

void init_lattice(int *particles, int particles_size, int *lat, int lat_size, 
        int species) {
    for (int i = 0; i < particles_size; i++){
        int pos = particles[2*i];
        int spec = particles[2*i+1];
        lat[lat_size*spec + pos] += 1;
    }
}

double single_species_hop_rate(double n, double b, double gamma) {
    return 1. + b/n;
}

double my_rand_double(pcg32_random_t* rngptr){
    return ldexp(pcg32_random_r(rngptr), -32); 
}

int periodic_bc(int pos, int lat_size){
    if (pos >= lat_size) {
        pos -= lat_size;
    } else if (pos < 0) {
        pos += lat_size;
    }
    
    return pos;
}

void update_particle_position_periodic_bc(pcg32_random_t* rngptr, int* particles, 
        int particles_size, int particle_index, int particle_pos, 
        int particle_spec, int* lat, int lat_size){

    int step = 2*pcg32_boundedrand_r(rngptr, 2) - 1;
    //int step = 1;
    int new_pos = periodic_bc(particle_pos + step, lat_size);
    
    particles[particle_index] = new_pos;
    lat[particle_pos + lat_size*particle_spec] -= 1;
    lat[new_pos + lat_size*particle_spec] += 1;
}



void init_particles(pcg32_random_t* rngptr, int *arr, int L, int N, int S) {
    for(int i=0; i < 2*N; i+=2){
        //arr[i] = pcg32_boundedrand_r(rngptr, L);
        arr[i] = L / 2;
        arr[i+1] = pcg32_boundedrand_r(rngptr, S);
    }
}

void sim(pcg32_random_t* rngptr, int* particles, int particles_size, int* lat, 
        int lat_size, double dt, double b, double gamma) {
    int rand_particle_index;
    int rand_pos;
    int rand_spec;
    int lat_occ;
    double hop_prob;
    double hop_event;
    
    for (int i = 0; i < particles_size; i++) {
        rand_particle_index = 2*pcg32_boundedrand_r(rngptr, particles_size);
        rand_pos = particles[rand_particle_index];
        rand_spec = particles[rand_particle_index + 1];
        lat_occ = lat[rand_pos + lat_size*rand_spec];

        hop_prob =  single_species_hop_rate((double)lat_occ, b, gamma)*dt;
        printf("%d, %lf, %lf, ", lat_occ, hop_prob, dt);     
        hop_event = my_rand_double(rngptr);
        printf("%lf\n", hop_event);

        if (hop_event < hop_prob) {
            update_particle_position_periodic_bc(rngptr, particles, 
                particles_size, rand_particle_index, rand_pos, rand_spec, lat, 
                lat_size);
        }
        
    }
}

void sim_test(pcg32_random_t* rngptr, int* particles, int particles_size, int* lat, 
        int lat_size, double dt, double b) {
        for (int i = 0; i < 2*particles_size; i++)
            particles[i] = 1;
    }

