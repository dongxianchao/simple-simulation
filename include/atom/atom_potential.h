#include "../structure/structure.h"
#include "../coefficient/coefficient.h"
#include <cstdlib>
#include <array>
#include <fstream>

class atom
{
    private:

    void transform_coordinate(cell &cel);

    double find_temperature(double& energy);

    double potentialenergy;
    double kineticenergy; 




    public:

    std::vector<std::vector<double>> basic_vectors; // it's used for describe the structure of cell
    std::vector<std::vector<double>> basic_vectors_inv;
    
    
    double box[6];
    int boxes_x, boxes_y, boxes_z; // numbers of the neighbour box in 3 dim
    int* boxes;
    int* neighbour_list;

    int total_num;
    int cell_num;
    double max_nl_dr_square = 0;
    std::vector<double> x, y, z, mass, vx, vy, vz, pe, fx, fy, fz, mass_inv, x0, y0, z0, rho;


    atom(std::vector<int> numcells, int numbers,
    std::vector<double> mass_0,
    std::vector<double> x_0,
    std::vector<double> y_0,
    std::vector<double> z_0, cell& cell_0); /*give the num of cells in 3
    dimensions and the numbers, placement of atoms in the each cell*/

    double kinetic_energy();// this is used for calculate the total energy
    double total_energy();
    double potential_energy();


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



void EAM_cubic_spline_interpolation_partial(int& num, double*& values, double& step);



class potential_energy
{
    private:

    friend void EAM_cubic_spline_interpolation_partial(int& num, double*& values, double& step);

    static double EAM_drho;
    static double EAM_dr;
    static int EAM_Nrho;
    static int EAM_Nr;
    static double EAM_drho_inv;
    static double EAM_dr_inv;
    static double* EAM_F_rho;
    static double* EAM_Z_r;
    static double* EAM_rho_r;
    static double* EAM_cubic_spline_coefficient_r_Zr;
    static double* EAM_cubic_spline_coefficient_rho_Frho;   
    static double* EAM_cubic_spline_coefficient_r_rhor;

    public:

    static void LJ_potential(atom & atom_0);    //it's classical method
    static void LJ_potential_nl(atom & atom_0);
    // it's the near boxes method

    
    static void get_the_EAM_data(std::ifstream& file);
    static void EAM_potential_cubicspline(atom& atom_0);  //the same as the above
    static void EAM_cubic_spline_interpolation();   // used for getting the parameter
    static void EAM_nl_linear_interpolation(atom& atom_0);  // used for update force and pe
    


    static double get_the_value_of_spline(double*& coefficient_num, double& x, double& step_inv);
    static double get_the_derivative_of_the_spline(double*& coefficient_num, double& x, double& step_inv);
};



