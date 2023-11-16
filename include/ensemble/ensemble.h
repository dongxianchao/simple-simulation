#include "../atom/atom_potential.h"
#include <random>
class ensemble
{
public:
static std::string ensemble_thermostat;
static double T0_Nf;
static double e_index_couple;
static double e_index_couple_root;//    above variables are used for lan 
static int ber_pressure_condition;  //this parameter is used for determine the type of ber NPT

static double pressure_hydro;
static double beta_hydro;   //parameters for condition 1

static double pressure_xx;
static double pressure_xy;
static double pressure_xz;
static double pressure_yx;
static double pressure_yy;
static double pressure_yz;
static double pressure_zx;
static double pressure_zy;
static double pressure_zz; 

static double beta_xx;
static double beta_xy;
static double beta_xz;
static double beta_yx;
static double beta_yy;
static double beta_yz;
static double beta_zx;
static double beta_zy;
static double beta_zz;


//parameters for condition 2 and 3



static std::default_random_engine engine;
static double target_T, T_coupling_coefficient, T_coupling_coefficient_inv, P_coupling_coefficient, P_coupling_coefficient_inv;
static std::normal_distribution<double>* gaussian;
static std::chi_squared_distribution<double>* chi_squared;

//this is used for NVT ensemble, usually operated after algorithms
static void NVT_ensemble(atom& atom_0);


//this is used for NPT ensemble, usually operated after algorithms
static void NPT_ensemble(atom& atom_0);
};
