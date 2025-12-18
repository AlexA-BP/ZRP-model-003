#include <stdio.h>
#include <stdlib.h>
//#include "simulation.h"
#include "pcg_variants.h"
#include "entropy.h"

typedef struct {
    const int N;
    const int L;
    const int species;
    const double b;
    const char* bc;
    const double dt;
    const int time;
} params_t;

typedef struct{
    int *arr;
    int arr_c; 
} arr_int_t;

typedef struct{
    double *arr;
    int arr_c;
} arr_double_t;



int main(){

    const int N = 1000;
    const int L = 300;
    const int species = 1;
    const double b = 3.5;
    const char* bc = "periodic";
    const double dt = 1./(1.+b);
    const int time = 1000;

    params_t params = {
        .N = N, 
        .L = L,
        .species = species,
        .b = b,
        .bc = bc,
        .dt = dt,
        .time = time 
    };

    // setting up rng
    pcg32_random_t rng; 
    uint64_t seeds[2];
    entropy_getbytes((void*)seeds, sizeof(seeds));
    pcg32_srandom_r(&rng, seeds[0], seeds[1]);

    FILE *file = fopen("./data/test.dat", "w");
    if (file == NULL){
        printf("Error opening file!\n");
        return 1;
    }

    int* particles = calloc(N * 2, sizeof(int));
    if (!particles) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    int* lat = calloc(params.L * params.species, sizeof(int));
    if (!lat) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    init_particles(&rng, particles, params.L, params.N, params.species);
    init_lattice(particles, params.N, lat, params.L, params.species);


    for (int i = 0; i < params.time; i++){
        sim(&rng, particles, N, lat, L, dt, b, 1.);

        fprintf(file, "t = %d: ", i);
        for (int k = 0; k < L; k++){
            fprintf(file, "%d ", lat[k]);    
        }
        fprintf(file, "\n");   
    }

    fclose(file);

    return 0;
}