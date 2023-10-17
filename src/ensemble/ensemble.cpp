#include "../../include/ensemble/ensemble.h"
#include <iostream>



double ensemble::T0_Nf;
double ensemble::e_index_couple;
double ensemble::e_index_couple_root;
std::string ensemble::ensemble_thermostat;
double ensemble::target_T, ensemble::coupling_coefficient, ensemble::coupling_coefficient_inv;
std::random_device rd;
std::default_random_engine ensemble::engine(rd()); 
std::normal_distribution<double>* ensemble::gaussian;
std::chi_squared_distribution<double>* ensemble::chi_squared;

void ensemble::NVT_ensemble(atom& atom_0)
{
    if (ensemble_thermostat == "ber")
    {
        double scalar = sqrt(1 + coupling_coefficient_inv * (target_T / atom_0.T - 1));
        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
            atom_0.velocity[i].x *= scalar;
            atom_0.velocity[i].y *= scalar;
            atom_0.velocity[i].z *= scalar;
        }
        atom_0.kinetic_energy();
        atom_0.find_temperature();
        atom_0.total_energy();  //update the information
    }
    if (ensemble_thermostat == "bdp")
    {
        double R = (*gaussian)(engine);
        double sum_R_square = (*chi_squared)(engine);
        double scalar = sqrt(e_index_couple + T0_Nf / atom_0.T * (1 - e_index_couple)
        * (R * R + sum_R_square) + 2 * e_index_couple_root * R * sqrt(T0_Nf / atom_0.T * (1 - e_index_couple)));
        //  get the scalar
        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
            atom_0.velocity[i].x *= scalar;
            atom_0.velocity[i].y *= scalar;
            atom_0.velocity[i].z *= scalar;
        }
        atom_0.kinetic_energy();
        atom_0.find_temperature();
        atom_0.total_energy();  //update the information
    }
}