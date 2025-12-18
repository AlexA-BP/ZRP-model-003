#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "hdf5.h"
#include "simulation.h"
#include "pcg_variants.h"
#include "entropy.h"

#define FILENAME "data/dset.h5"

void H5_add_int_attributes(hid_t h5id, hid_t attr_type, const char* attr_name, 
    const int attr_val);

void progress(long double part, long double tot); 

int main(){ 
    const int TIME = 5 * 1000 * 5; // t/tau * tau * 1/dt
    const int N = 2*1230;
    const int L = 1230;
    const int SPECIES = 1;
    const double gamma = 1.;
    const double b = 4.;
    const double dt = 1./(b);
    const int chunk_size = 10 * 1000;

    
    // setting up rng
    pcg32_random_t rng; 
    uint64_t seeds[2];
    entropy_getbytes((void*)seeds, sizeof(seeds));
    pcg32_srandom_r(&rng, seeds[0], seeds[1]);


    // setting up HDF5
    hid_t file_id, dataspace_id, plist_id, dataset_id;
    herr_t status __attribute__((unused));
    hsize_t dims[3] = {TIME, N, 2};
    hsize_t chunk_dims[3] = {chunk_size, N, 2};

    
    file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating HDF5 file\n");
        return 1;
    }

    // Create dataspace (unlimited size for first dimension)
    dataspace_id = H5Screate_simple(3, dims, NULL);
    if (dataspace_id < 0) {
        fprintf(stderr, "Error creating dataspace\n");
        H5Fclose(file_id);
        return 1;
    }

    // Create dataset with chunking for efficient streaming
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    if (plist_id == H5I_INVALID_HID){
        fprintf(stderr, "Error creating chunking\n");
        H5Fclose(file_id);
        H5Sclose(dataspace_id);
    }

    H5Pset_chunk(plist_id, 3, chunk_dims);
    H5Pset_deflate(plist_id, 6);  // compression

    dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_INT, dataspace_id,
                                  H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID){
        fprintf(stderr, "Error creating dataset\n");
        H5Pclose(plist_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
    }

    // Allocate memory for one array (M × L)
    int* particles = calloc(N * 2, sizeof(int));
    if (!particles) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    int* lat = calloc(L * SPECIES, sizeof(int));
    if (!lat) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    init_particles(&rng, particles, L, N, SPECIES);
    init_lattice(particles, N, lat, L, SPECIES);

    // === Add attributes
    // H5_add_int_attributes(dataset_id, H5T_NATIVE_INT, "N", N);
    // H5_add_int_attributes(dataset_id, H5T_NATIVE_INT, "M", M);
    // H5_add_int_attributes(dataset_id, H5T_NATIVE_INT, "L", L);

    // Write arrays in chunks
    for (int i = 0; i < TIME; i += chunk_size){
        int current_chunk_size = 
            (i + chunk_size > TIME) ? (TIME - i) : chunk_size;

        // Create hyperslab in dataset
        hsize_t offset[3] = {i, 0, 0};
        hsize_t count[3] = {current_chunk_size, N, 2};

        // Select hyperslab in dataset
        status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, 
            count, NULL);

        // Create memory dataspace for the chunk
        hid_t mem_dataspace_id = H5Screate_simple(3, count, NULL);
        
        // Fill the chunk data (one array at a time)
        int* chunk_data = malloc(current_chunk_size * N * 2 * sizeof(int));
        if (chunk_data == NULL){
            fprintf(stderr, "Memory allocation failed for the chunk\n");
            H5Sclose(mem_dataspace_id);
            H5Dclose(dataset_id);
            H5Pclose(plist_id);
            H5Sclose(dataspace_id);
            H5Fclose(file_id);
            free(particles);
            free(lat);
            return 1;
        }

        // Fill chunk_data with simulation results:
        for (int j = 0; j < current_chunk_size; j++){
            sim(&rng, particles, N, lat, L, dt, b, gamma);

            for (int k = 0; k < N * 2; k++){
                chunk_data[j*N*2 + k] = particles[k];
                // printf("\n k=%d, particles[k]=%d\n", k, particles[k]);
                //chunk_data[j*N*2 + k] = 1;
            }
        }
        

        // Write chunk to HDF5
        status = H5Dwrite(dataset_id, H5T_NATIVE_INT, mem_dataspace_id, 
            dataspace_id, H5P_DEFAULT, chunk_data);
        

        if (status == H5I_INVALID_HID){
            fprintf(stderr, "Error writing chunk %d\n", i);
        } else {
            progress((long double)i, (long double)TIME);
        }

        // cleanup 
        free(chunk_data);
        H5Sclose(mem_dataspace_id);

    } 
    printf("\rProgress (100.00%%)");
    fflush(stdout);
    

    // cleanup
    H5Dclose(dataset_id);
    H5Pclose(plist_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);


    free(particles);
    free(lat);

    printf("\nExecution complete\n");

    return 0;
}

void progress(long double part, long double tot) {
    long double perc = 100 * part / tot;
    printf("\rProgress (%.2Lf%%)", perc);
    fflush(stdout); 
}


void H5_add_int_attributes(hid_t h5id, hid_t attr_type, const char* attr_name, 
    int attr_val){
    hsize_t attr_dims = 1;
    hid_t attr_type_cp = H5Tcopy(attr_type);
    hid_t attr_space = H5Screate_simple(1, &attr_dims, NULL);
    hid_t attr_id = H5Acreate2(h5id, attr_name, attr_type_cp, attr_space,
                               H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id != H5I_INVALID_HID) {
        H5Awrite(attr_id, attr_type, &attr_val);
        H5Aclose(attr_id);
    }
    H5Sclose(attr_space);
    H5Tclose(attr_type_cp);
}