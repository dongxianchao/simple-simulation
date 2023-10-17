#pragma once
#include "../structure/structure.h"
#include "../coefficient/coefficient.h"
#include <cstdlib>
#include <array>
#include <fstream>

struct three_dim_vector
{
    double x = 0, y = 0, z = 0;
};


class atom
{
    private:

    void transform_coordinate(cell &cel);






    public:

    std::vector<std::vector<double>> basic_vectors; // it's used for describe the structure of cell
    std::vector<std::vector<double>> basic_vectors_inv;

    double potentialenergy; //this is the potentialenergy
    double kineticenergy; //    this is the kineticenergy
    double totalenergy;//   this is the total energy
    double T;   //this is the temprature

    

    double box[6];
    int boxes_x, boxes_y, boxes_z; // numbers of the neighbour box in 3 dim
    int* boxes;
    int* neighbour_list;

    int total_num;
    int cell_num;
    double max_nl_dr_square = 0;
    three_dim_vector *coordinate, *velocity, *force, *coordinate_init;
    std::vector<double> mass, pe, mass_inv, rho;
    std::vector<int> type, count;
    bool if_othogonal;
    


    atom(std::vector<int> numcells, int numbers,
    std::vector<double> mass_0,
    std::vector<double> x_0,
    std::vector<double> y_0,
    std::vector<double> z_0, cell& cell_0); /*give the num of cells in 3
    dimensions and the numbers, placement of atoms in the each cell*/

    double kinetic_energy();    // this is used for calculate the total energy
    double total_energy();
    double potential_energy();
    void update_energy();   //this method is used for updating the energy information
    double find_temperature();

    void initialise_boxes();//  used to initialise boxes to empty
    void update_boxes();// used to update the boxes

    void initialise_neighbour_list();// used to initialise_neighbour_list
    void update_neighbour_list();   //  used to update




    void initialise_velocity(double T0);

    std::vector<std::vector<double>> print_status();  //print the output

    void initialise_boxes_xyz_maxnums(cell& cell_0);//  used for neighbour algorithm

    int *get_place(double& x, double& y, double& z, atom& atom_0);//  used for get index for boxes

    void insert_in_boxes(int*& place, int index);
    void insert_in_boxes(int (&place)[3], int index);

};

void compute_1dv_times_3dm(double* arr, std::vector<std::vector<double>> &mtx, double (&arr_out)[3]);







class potential_energy
{
    private:

    static double EAM_drho;//   the interval of rho
    static double EAM_dr;// the interval of r
    static int EAM_Nrho;//  the total number of points about rho
    static int EAM_Nr;// the total number of points about r
    static double EAM_drho_inv;//   the inverse of the interval of rho
    static double EAM_dr_inv;// the inverse of the interval of r
    static double* EAM_F_rho;// the list of F about rho
    static double* EAM_Z_r;//   the list of  Z of r
    static double* EAM_rho_r;   // the list of rho of r
    static double* EAM_cubic_spline_coefficient_r_Zr;// the parameters of cubic spline
    static double* EAM_cubic_spline_coefficient_rho_Frho;   
    static double* EAM_cubic_spline_coefficient_r_rhor;
    static double** EAM_alloy_F_rho;
    static double** EAM_alloy_rfhi_r;
    static double** EAM_alloy_rho_r;
    static double** EAM_alloy_cubic_spline_coefficient_rfhi_r;
    static double** EAM_alloy_cubic_spline_coefficient_rho_r;
    static double** EAM_alloy_cubic_spline_coefficient_F_rho;// the parameters of cubic spline
    static int total_pair_nums;
    static int atom_types_nums;



    public:

    static void LJ_potential(atom & atom_0);    //it's classical method
    static void LJ_potential_nl(atom & atom_0);
    // it's the near boxes method

    
    static void get_the_EAM_data(std::ifstream& file);  //to generate the data
    static void EAM_potential_cubicspline(atom& atom_0);  //update the force
    static void EAM_cubic_spline_interpolation();   // used for getting the parameter
    static void EAM_nl_linear_interpolation(atom& atom_0);  // used for update force and pe

    


    static void get_the_EAM_alloy_data(std::ifstream& file, atom& atom_0, bool random=false, int type=0, double percentage=0);// used for generate the data
    static void EAM_alloy_cubic_spline_interpolation();//   used for getting the parameter
    static void EAM_alloy_potential_cubicspline(atom& atom_0);// update the force
    static void EAM_alloy_nl_cubic_spline(atom& atom_0);    //used for update force and pe



    static double get_the_value_of_spline(double*& coefficient_num, double& x, double& step_inv);
    static double get_the_derivative_of_the_spline(double*& coefficient_num, double& x, double& step_inv);
};



