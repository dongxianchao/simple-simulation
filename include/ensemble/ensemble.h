#include "../atom/atom_potential.h"
#include <random>
class ensemble
{
public:
static std::string ensemble_thermostat;
static double T0_Nf;
static double e_index_couple;
static double e_index_couple_root;//    above variables are used for lan 
static std::default_random_engine engine;
static double target_T, coupling_coefficient, coupling_coefficient_inv;
static std::normal_distribution<double>* gaussian;
static std::chi_squared_distribution<double>* chi_squared;
static void NVT_ensemble(atom& atom_0);
};
