#include "../../include/ensemble/ensemble.h"
#include <iostream>



double ensemble::T0_Nf;
double ensemble::e_index_couple;
double ensemble::e_index_couple_root;
int ensemble::ber_pressure_condition = 0;

double ensemble::beta_hydro = 0;
double ensemble::pressure_hydro = 0;

double ensemble::beta_xx = 0;
double ensemble::beta_xy = 0;
double ensemble::beta_xz = 0;
double ensemble::beta_yx = 0;
double ensemble::beta_yy = 0;
double ensemble::beta_yz = 0;
double ensemble::beta_zx = 0;
double ensemble::beta_zy = 0;
double ensemble::beta_zz = 0;
double ensemble::pressure_xx = 0;
double ensemble::pressure_xy = 0;
double ensemble::pressure_xz = 0;
double ensemble::pressure_yx = 0;
double ensemble::pressure_yy = 0;
double ensemble::pressure_yz = 0;
double ensemble::pressure_zx = 0;
double ensemble::pressure_zy = 0;
double ensemble::pressure_zz = 0;

std::string ensemble::ensemble_thermostat;
double ensemble::target_T, ensemble::T_coupling_coefficient, ensemble::T_coupling_coefficient_inv, ensemble::P_coupling_coefficient, ensemble::P_coupling_coefficient_inv;
std::random_device rd;
std::default_random_engine ensemble::engine(rd()); 
std::normal_distribution<double>* ensemble::gaussian;
std::chi_squared_distribution<double>* ensemble::chi_squared;

void ensemble::NVT_ensemble(atom& atom_0)
{
    if (ensemble_thermostat == "ber")
    {
        double scalar = sqrt(1 + T_coupling_coefficient_inv * (target_T / atom_0.T - 1));
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

void ensemble::NPT_ensemble(atom& atom_0)
{
    if (ensemble_thermostat == "ber")
    {
        NVT_ensemble(atom_0);   
        atom_0.calculate_virial_pressure();
        if (ber_pressure_condition == 1)
        {
            double p_real = (atom_0.virial_press_tensor_tot[0] + atom_0.virial_press_tensor_tot[3] + atom_0.virial_press_tensor_tot[6]) / 3;
            double miu_hydro = 1 - beta_hydro * P_coupling_coefficient_inv / 3 * (pressure_hydro - p_real);
            atom_0.rescale_the_position(miu_hydro, miu_hydro, miu_hydro);
        }

        if (ber_pressure_condition == 2)
        {
            double miu_xx = 1 - beta_xx * P_coupling_coefficient_inv / 3 * (pressure_xx - atom_0.virial_press_tensor_tot[0]);
            double miu_yy = 1 - beta_yy * P_coupling_coefficient_inv / 3 * (pressure_yy - atom_0.virial_press_tensor_tot[3]);
            double miu_zz = 1 - beta_zz * P_coupling_coefficient_inv / 3 * (pressure_zz - atom_0.virial_press_tensor_tot[6]);
            atom_0.rescale_the_position(miu_xx, miu_yy, miu_zz);
        }
        if (ber_pressure_condition == 3)
        {
            double miu[9]{};
            miu[0] = 1 - beta_xx * P_coupling_coefficient_inv / 3 * (pressure_xx - atom_0.virial_press_tensor_tot[0]);
            miu[1] = 1 - beta_xy * P_coupling_coefficient_inv / 3 * (pressure_xy - atom_0.virial_press_tensor_tot[1]);
            miu[2] = 1 - beta_xz * P_coupling_coefficient_inv / 3 * (pressure_xz - atom_0.virial_press_tensor_tot[2]);
            miu[3] = 1 - beta_yx * P_coupling_coefficient_inv / 3 * (pressure_yx - atom_0.virial_press_tensor_tot[3]);
            miu[4] = 1 - beta_yy * P_coupling_coefficient_inv / 3 * (pressure_yy - atom_0.virial_press_tensor_tot[4]);
            miu[5] = 1 - beta_yz * P_coupling_coefficient_inv / 3 * (pressure_yz - atom_0.virial_press_tensor_tot[5]);
            miu[6] = 1 - beta_zx * P_coupling_coefficient_inv / 3 * (pressure_zx - atom_0.virial_press_tensor_tot[6]);
            miu[7] = 1 - beta_zy * P_coupling_coefficient_inv / 3 * (pressure_zy - atom_0.virial_press_tensor_tot[7]);
            miu[8] = 1 - beta_zz * P_coupling_coefficient_inv / 3 * (pressure_zz - atom_0.virial_press_tensor_tot[8]);
            atom_0.rescale_the_position(miu);
        }
    }
}