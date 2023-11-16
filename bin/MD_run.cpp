#include <fstream>
#include <iostream>
#include "../include/algorithm/algorithm.h"
#include "../include/ensemble/ensemble.h"

template <typename T>
void get_vector(std::vector<T>& v, std::ifstream & file, int n)
{
    v.resize(n);
    for (int i = 0; i < n; ++i)
    {
        file >> v[i];
    }
}


void output_xyz(atom*& atom_0, std::fstream& file, int& step)
{
    file << atom_0->total_num << "\n"; // out put the number
    file << "step" << "\t" << step << "\n";     //out put the steps
    for (int i = 0; i < atom_0->total_num; ++i)
    {
        double x = atom_0->coordinate[i].x;
        double y = atom_0->coordinate[i].y;
        double z = atom_0->coordinate[i].z;
        if (!atom_0->if_othogonal)
        {
            double old_r[3] = {x, y, z};
            double new_r[3];
            compute_1dv_times_3dm(old_r, atom_0->basic_vectors, new_r);
            x = new_r[0];
            y = new_r[1];
            z = new_r[2];
        }
        file << atom_0->type[i] << "\t" << x << "\t" << y << "\t" << z << "\n"; // out put the coordinates
    }
}

cell* handle_cell(std::ifstream& file)
{
    std::vector<double> cell_x(3), cell_y(3), cell_z(3), lattice(3);
    get_vector(cell_x, file, 3);
    get_vector(cell_y, file, 3);
    get_vector(cell_z, file, 3);
    get_vector(lattice, file, 3);
    std::vector<std::vector<double>> cell_basic = {cell_x, cell_y, cell_z};
    cell* cell_0 = new cell(lattice[0], lattice[1], lattice[2], cell_basic);
    return cell_0;
}

atom* handle_atom(std::ifstream& file, cell* cell_0)
{
    int total;
    file >> total;
    std::vector<int> numcells(total);
    std::vector<double> x(total), y(total), z(total), mass(total);
    get_vector(numcells, file, 3);
    get_vector(x, file, total);
    get_vector(y, file, total);
    get_vector(z, file, total);
    get_vector(mass, file, total);
    atom* atom_0 = new atom(numcells, total, mass, x, y, z, *cell_0);
    return atom_0;
}

atom* handle_atom(std::ifstream& file, cell* cell_0, atom*& atom_1)//   used for test
{
    int total;
    file >> total;
    std::vector<int> numcells(total);
    std::vector<double> x(total), y(total), z(total), mass(total);
    get_vector(numcells, file, 3);
    get_vector(x, file, total);
    get_vector(y, file, total);
    get_vector(z, file, total);
    get_vector(mass, file, total);
    atom* atom_0 = new atom(numcells, total, mass, x, y, z, *cell_0);
    atom_1 = new atom(numcells, total, mass, x, y, z, *cell_0);
    return atom_0;
}

void handle_ensemble(std::ifstream& file, std::string& ensemble_type, atom*& atom_0)
{
    file >> ensemble_type;
    if (ensemble_type == "nvt")
    {
        file >> ensemble::ensemble_thermostat >> ensemble::target_T;
        if (ensemble::ensemble_thermostat == "ber")
        {
            file >> ensemble::T_coupling_coefficient;  //read the nvt ensemble 
            ensemble::T_coupling_coefficient_inv = 1 / ensemble::T_coupling_coefficient;
        }
        if (ensemble::ensemble_thermostat == "bdp")
        {
            file >> ensemble::T_coupling_coefficient;  //read the nvt ensemble 
            ensemble::T_coupling_coefficient_inv = 1 / ensemble::T_coupling_coefficient;
            double kinetic_energy;
            kinetic_energy = 1.5 * coefficient::k_B * atom_0->total_num * ensemble::target_T;
            double Nf = 3 * atom_0->total_num;
            ensemble::T0_Nf = ensemble::target_T / Nf;
            ensemble::e_index_couple = exp(-ensemble::T_coupling_coefficient_inv);
            ensemble::e_index_couple_root = exp(-0.5*ensemble::T_coupling_coefficient_inv);
            ensemble::gaussian = new std::normal_distribution<double>(0, 1);
            ensemble::chi_squared = new std::chi_squared_distribution<double>(3 * atom_0->total_num - 1);
        }
    }
    if (ensemble_type == "npt")
    {
        atom_0->virial_switch = true;
        file >> ensemble::ensemble_thermostat >> ensemble::target_T;
        if (ensemble::ensemble_thermostat == "ber")
        {
            file >> ensemble::T_coupling_coefficient;
            ensemble::T_coupling_coefficient_inv = 1 / ensemble::T_coupling_coefficient;
            file >> ensemble::ber_pressure_condition;
            switch (ensemble::ber_pressure_condition)
            {
                case 1:
                file >> ensemble::pressure_hydro >> ensemble::beta_hydro >> ensemble::P_coupling_coefficient;
                ensemble::P_coupling_coefficient_inv = 1 / ensemble::P_coupling_coefficient;
                ensemble::pressure_hydro /= coefficient::unit_p_bar;
                ensemble::beta_hydro *= coefficient::unit_p_bar;
                break;

                case 2:
                file >> ensemble::pressure_xx >> ensemble::pressure_yy >> ensemble::pressure_zz >> ensemble::beta_xx >> ensemble::beta_yy >> ensemble::beta_zz
                >> ensemble::P_coupling_coefficient;
                ensemble::P_coupling_coefficient_inv = 1 / ensemble::P_coupling_coefficient;
                break;

                case 3:
                file >> ensemble::pressure_xx >> ensemble::pressure_xy >> ensemble::pressure_xz
                >> ensemble::pressure_yx >> ensemble::pressure_yy >> ensemble::pressure_yz
                >> ensemble::pressure_zx >> ensemble::pressure_zy >> ensemble::pressure_zz;
                file >> ensemble::beta_xx >> ensemble::beta_xy >> ensemble::beta_xz
                >> ensemble::beta_yz >> ensemble::beta_yy >> ensemble::beta_yz
                >> ensemble::beta_zx >> ensemble::beta_zy >> ensemble::beta_zz;
                file >> ensemble::P_coupling_coefficient;
                ensemble::P_coupling_coefficient_inv = 1 / ensemble::P_coupling_coefficient;
                break;
            }
        }
    }

}

void apply_ensemble(atom* atom_0, std::string& ensemble)
{
    if (ensemble == "nvt")
    ensemble::NVT_ensemble(*atom_0);
    if (ensemble == "npt")
    ensemble::NPT_ensemble(*atom_0);
}


void handle_parameter(std::ifstream& file, double& T, double& dt, double& step)
{
    std::string rubbish;    //used for skip the first word
    file >> rubbish >> T >> rubbish >> step >> rubbish >> dt;   //get the parameter
}

void handle_LJ_potential(std::ifstream& file)
{
    std::string rubbish;
    double epsilon, sigma, cutoff;
    file >> rubbish >> epsilon >> rubbish >> sigma >> rubbish >> cutoff;
    coefficient::get_LJ_coefficient(epsilon, sigma, cutoff);
}

void run_by_LJ_potential(atom*& atom_0, cell*& cell_0, double step, double dt, std::string ensemble)
{
    double dt_inner = dt / coefficient::unit_t_ps;
    atom_0 -> initialise_boxes_xyz_maxnums(*cell_0);
    int step_to_record = (0.01 * step > 1) ? 0.01 * step : 1;
    int step_to_print = (0.1 * step > 1) ? 0.1 * step : 1;
    if (atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4)
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_LJ(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);

            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
        }
    }

    else
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_LJ_nl(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);
            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
    }
    }   
}
void handle_EAM_potential()
{

    std::ifstream file("EAM.txt");
    if(!file)
    {
        std::cerr << "Error: unable to open file" << std::endl;
    }

    potential_energy::get_the_EAM_data(file);


    
}

void run_by_EAM_potential(atom*& atom_0, cell*& cell_0, double step, double dt, std::string ensemble)
{
    double dt_inner = dt / coefficient::unit_t_ps;
    atom_0 -> initialise_boxes_xyz_maxnums(*cell_0);
    potential_energy::EAM_cubic_spline_interpolation();//    used for get parameter
    int step_to_record = (0.01 * step > 1) ? 0.01 * step : 1;
    int step_to_print = (0.1 * step > 1) ? 0.1 * step : 1;
    if (atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4)//atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_EAM(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);

            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                static std::fstream xyzfile("../output/coordinates.xyz", std::ios::out);
                outfile << atom_0->T << "\t" << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
        }
    }

    else
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_EAM_nl(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);
            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                static std::fstream xyzfile("../output/coordinates.xyz", std::ios::out);
                outfile << atom_0->T << "\t" << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
                output_xyz(atom_0, xyzfile, i);
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
    }
    }   
}


void handle_EAM_alloy_potential(atom& atom_0, std::ifstream& file_0)
{
    bool random = false;
    double percentage;
    int type;
    std::string command;
    file_0 >> command;
    if (command == "random")
    {
        random = true;
        file_0 >> type >> percentage;
    }

    std::ifstream file("EAM_alloy.txt");
    if(!file)
    {
        std::cerr << "Error: unable to open file" << std::endl;
    }

    potential_energy::get_the_EAM_alloy_data(file, atom_0, random, type, percentage);
    
}

void run_by_EAM_alloy_potential(atom*& atom_0, cell*& cell_0, double step, double dt, std::string ensemble)
{
    double dt_inner = dt / coefficient::unit_t_ps;
    atom_0 -> initialise_boxes_xyz_maxnums(*cell_0);
    potential_energy::EAM_alloy_cubic_spline_interpolation();//    used for get parameter
    int step_to_record = (0.001 * step > 1) ? 0.001 * step : 1;
    int step_to_print = (0.1 * step > 1) ? 0.1 * step : 1;
    if (atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4)//atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_EAM_alloy(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);
            if (i % step_to_record == 0)//i % step_to_record == 0
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                static std::fstream xyzfile("../output/coordinates.xyz", std::ios::out);
                outfile << atom_0->T << "\t" << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
                output_xyz(atom_0, xyzfile, i);
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
        }
    }

    else
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_EAM_alloy_nl(dt_inner, *atom_0);
            atom_0->update_energy();
            apply_ensemble(atom_0, ensemble);
            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                static std::fstream xyzfile("../output/coordinates.xyz", std::ios::out);
                outfile << atom_0->T << "\t" << atom_0 -> potentialenergy << "\t" << atom_0 -> kineticenergy << "\t" 
                << atom_0 -> totalenergy << "\t" << (i * dt * 0.001) << "\n";
                output_xyz(atom_0, xyzfile, i);
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
    }
    }   
}

void test_for_EAM_alloy(atom*& atom_0, atom*& atom_1, double dt, cell*& cell_0)
{
    for (int i = 0; i < atom_0->total_num; ++i)
    {
        atom_1->velocity[i].x = atom_0->velocity[i].x;
        atom_1->velocity[i].y = atom_0->velocity[i].y;
        atom_1->velocity[i].z = atom_0->velocity[i].z;
        atom_1->mass[i] = atom_0->mass[i];
        atom_1->mass_inv[i] = atom_0->mass_inv[i];
        atom_1->type[i] = atom_0->type[i];
    }
    double dt_inner = dt / coefficient::unit_t_ps;
    atom_0 -> initialise_boxes_xyz_maxnums(*cell_0);
    potential_energy::EAM_alloy_cubic_spline_interpolation();//    used for get parameter

    for (int n = 0; n < 1; ++n)
    {
        algorithm::velvet_EAM_alloy_nl(dt_inner, *atom_0);
        algorithm::velvet_EAM_alloy(dt_inner, *atom_1);
        std::cout << atom_0->potential_energy() << "\t" << atom_1->potential_energy() << "\n";
    }

    std::fstream file("test.txt", std::ios::out);

    for (int i = 0; i < atom_0->total_num; ++i)
    {
        file << atom_0->type[i] << "\t" << atom_1->type[i] << "\n";
    }

}



int main()
{
    auto start_time = std::chrono::high_resolution_clock::now();


    std::ifstream file("INCAR.txt");

    if(!file)
    {
        std::cerr << "Error: unable to open file" << std::endl;
    }

    std::string line;

    cell* cell_0;

    atom* atom_0, *atom_1;

    double T, dt, step;

    std::string ensemble_type, ensemble_thermostat;


    while (std::getline(file, line))
    {
        if (line == "cell")
            cell_0 = handle_cell(file);      // read the cell
        if (line == "atom")
            atom_0 = handle_atom(file, cell_0, atom_1);     // read the atom
        if (line == "ensemble")
        {
            handle_ensemble(file, ensemble_type, atom_0);
        }
        if (line == "parameter")
        {
            handle_parameter(file, T, dt, step);    // read the parameter
            atom_0 -> initialise_velocity(T);
            break;
        }
    }


    while (std::getline(file, line))
    {
        if(line == "LJ potential")
        {
            handle_LJ_potential(file);
            run_by_LJ_potential(atom_0, cell_0, step, dt, ensemble_type);
        }

        if(line == "EAM potential")
        {
            handle_EAM_potential();
            run_by_EAM_potential(atom_0, cell_0, step, dt, ensemble_type);
        }

        if(line == "EAM_alloy potential")
        {
            handle_EAM_alloy_potential(*atom_0, file);
            atom_0->initialise_velocity(T);
            run_by_EAM_alloy_potential(atom_0, cell_0, step, dt, ensemble_type);
        }

    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Time taken: " << elapsed.count() << " ms" << std::endl;
    
    return 0;
}




