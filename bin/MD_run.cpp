#include <fstream>
#include <iostream>
#include "../include/algorithm/algorithm.h"
#include <chrono>
template <typename T>
void get_vector(std::vector<T>& v, std::ifstream & file, int n)
{
    v.resize(n);
    for (int i = 0; i < n; ++i)
    {
        file >> v[i];
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

void run_by_LJ_potential(atom*& atom_0, cell*& cell_0, double step, double dt)
{
    double dt_inner = dt / coefficient::unit_t_fs;
    atom_0 -> initialise_boxes_xyz_maxnums(*cell_0);
    int step_to_record = (0.01 * step > 1) ? 0.01 * step : 1;
    int step_to_print = (0.1 * step > 1) ? 0.1 * step : 1;
    if (atom_0 -> boxes_x < 4 || atom_0 -> boxes_y < 4 || atom_0 -> boxes_z < 4)
    {
        for (int i = 0; i < step; ++i)
        {
            algorithm::velvet_LJ(dt_inner, *atom_0);

            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potential_energy() << "\t" << atom_0 -> kinetic_energy() << "\t" 
                << atom_0 -> total_energy() << "\t" << (i * dt * 0.001) << "\n";
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
            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potential_energy() << "\t" << atom_0 -> kinetic_energy() << "\t" 
                << atom_0 -> total_energy() << "\t" << (i * dt * 0.001) << "\n";
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

void run_by_EAM_potential(atom*& atom_0, cell*& cell_0, double step, double dt)
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

            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potential_energy() << "\t" << atom_0 -> kinetic_energy() << "\t" 
                << atom_0 -> total_energy() << "\t" << (i * dt * 0.001) << "\n";
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
            if (i % step_to_record == 0)
            {
                static std::fstream outfile("../output/outfile.txt", std::ios::out);
                outfile << atom_0 -> potential_energy() << "\t" << atom_0 -> kinetic_energy() << "\t" 
                << atom_0 -> total_energy() << "\t" << (i * dt * 0.001) << "\n";
            }   // out put the file 

            if (i % step_to_print == 0)
            {
                std::cout << "still running" << std::endl;
            }
    }
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

    atom* atom_0;

    double T, dt, step;


    while (std::getline(file, line))
    {
        if (line == "cell")
            cell_0 = handle_cell(file);      // read the cell
        if (line == "atom")
            atom_0 = handle_atom(file, cell_0);     // read the atom
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
            run_by_LJ_potential(atom_0, cell_0, step, dt);
        }

        if(line == "EAM potential")
        {
            handle_EAM_potential();
            run_by_EAM_potential(atom_0, cell_0, step, dt);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Time taken: " << elapsed.count() << " ms" << std::endl;
    
    return 0;
}