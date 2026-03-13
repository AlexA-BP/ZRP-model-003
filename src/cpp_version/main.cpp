#include<iostream>

#include "H5Cpp.h"
#include "pcg_variants.h"
#include "entropy.h"

int main(int argc, char *argv[])
{
    std::cout << "You have entered " << argc
              << " arguments: " << std::endl;

    int i = 0;
    while (i < argc) 
    {
        std::cout << "Argument " << i + 1
                  << ": " << argv[i]
                  << std::endl;
        i++;
    }    

    return 0;
}