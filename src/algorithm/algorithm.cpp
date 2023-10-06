#include "../../include/algorithm/algorithm.h"
#include <iostream>
#include <algorithm>



//do some auxilary work
inline void mirror_constrain_once_copy(double& len,double& half_len, double &xij)
{
    if (xij > half_len)
    xij -= len;
    else if (xij < -half_len)
    xij += len;
}

void mirror_constrain_copy(double* b, double &xij, double &yij, double &zij)
// to apply the mirror constrain
{
    mirror_constrain_once_copy(b[0], b[3], xij);
    mirror_constrain_once_copy(b[1], b[4], yij);
    mirror_constrain_once_copy(b[2], b[5], zij);
}



inline void mirror_constrain_for_placement(atom &atom_0, double &x, double &y, double &z)
// to apply the mirror constrain
{
    while (x >= atom_0.box[0]) x -= atom_0.box[0];
    while (x < 0) x += atom_0.box[0];
    
    while (y >= atom_0.box[1]) y -= atom_0.box[1];
    while (y < 0) y += atom_0.box[1];

    while (z >= atom_0.box[2]) z -= atom_0.box[2];
    while (z < 0) z += atom_0.box[2];
}


// update the velocity

inline void algorithm::velvet_for_update_v(atom& atom_0, double& t)
{   
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i) {
        atom_0.velocity[i].x += 0.5 * (atom_0.force[i].x * atom_0.mass_inv[i]) * t; 
        atom_0.velocity[i].y += 0.5 * (atom_0.force[i].y * atom_0.mass_inv[i]) * t; 
        atom_0.velocity[i].z += 0.5 * (atom_0.force[i].z * atom_0.mass_inv[i]) * t; 
    }
}




//update the placement
inline void algorithm::velvet_for_update_place(double & dt, atom& atom_0)
{
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        double displacement_old[3];
        displacement_old[0] = atom_0.velocity[i].x * dt;
        displacement_old[1] = atom_0.velocity[i].y * dt;
        displacement_old[2] = atom_0.velocity[i].z * dt;       

        if (!atom_0.if_othogonal)
        {
            double displacement_new[3];
            compute_1dv_times_3dm(displacement_old, atom_0.basic_vectors_inv, displacement_new);
            displacement_old[0] = displacement_new[0];
            displacement_old[1] = displacement_new[1];
            displacement_old[2] = displacement_new[2];
        }

        atom_0.coordinate[i].x += displacement_old[0];
        atom_0.coordinate[i].y += displacement_old[1];
        atom_0.coordinate[i].z += displacement_old[2];
        mirror_constrain_for_placement(atom_0, atom_0.coordinate[i].x, atom_0.coordinate[i].y, atom_0.coordinate[i].z);

        if (std::isnan(atom_0.coordinate[i].x) || std::isnan(atom_0.coordinate[i].y) || std::isnan(atom_0.coordinate[i].z))
        throw std::out_of_range("the coordinate is wrong");
    }
}



//update the placement with neighbour list
inline void algorithm::velvet_for_update_place_with_nl(double & dt, atom& atom_0)
{
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        double displacement_old[3];
        displacement_old[0] = atom_0.velocity[i].x * dt;
        displacement_old[1] = atom_0.velocity[i].y * dt;
        displacement_old[2] = atom_0.velocity[i].z * dt;  
             
    
        if (!atom_0.if_othogonal)
        {
            double displacement_new[3];
            compute_1dv_times_3dm(displacement_old, atom_0.basic_vectors_inv, displacement_new);
            displacement_old[0] = displacement_new[0];
            displacement_old[1] = displacement_new[1];
            displacement_old[2] = displacement_new[2];
        }



        atom_0.coordinate[i].x += displacement_old[0];
        atom_0.coordinate[i].y += displacement_old[1];
        atom_0.coordinate[i].z += displacement_old[2];
        mirror_constrain_for_placement(atom_0, atom_0.coordinate[i].x, atom_0.coordinate[i].y, atom_0.coordinate[i].z);

        double dx, dy, dz;

        dx = atom_0.coordinate[i].x - atom_0.coordinate_init[i].x;
        dy = atom_0.coordinate[i].y - atom_0.coordinate_init[i].y;
        dz = atom_0.coordinate[i].z - atom_0.coordinate_init[i].z;

        mirror_constrain_copy(atom_0.box, dx, dy, dz);// prevent the atom pass through

        double dr_square = dx * dx + dy * dy + dz * dz;
        if (dr_square > atom_0.max_nl_dr_square)
        atom_0.max_nl_dr_square = dr_square;   

        
    }
}




inline void algorithm::velvet_for_update_v_first_step(atom& atom_0, double& t)
{
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i) {
        atom_0.velocity[i].x += 0.5 * (atom_0.force[i].x * atom_0.mass_inv[i]) * t; 
        atom_0.velocity[i].y += 0.5 * (atom_0.force[i].y * atom_0.mass_inv[i]) * t; 
        atom_0.velocity[i].z += 0.5 * (atom_0.force[i].z * atom_0.mass_inv[i]) * t;


        double displacement_old[3];
        displacement_old[0] = atom_0.velocity[i].x * t;
        displacement_old[1] = atom_0.velocity[i].y * t;
        displacement_old[2] = atom_0.velocity[i].z * t;  


             
    
        if (!atom_0.if_othogonal)
        {
            double displacement_new[3];
            compute_1dv_times_3dm(displacement_old, atom_0.basic_vectors_inv, displacement_new);
            displacement_old[0] = displacement_new[0];
            displacement_old[1] = displacement_new[1];
            displacement_old[2] = displacement_new[2];
        }




        atom_0.coordinate[i].x += displacement_old[0];
        atom_0.coordinate[i].y += displacement_old[1];
        atom_0.coordinate[i].z += displacement_old[2];
        mirror_constrain_for_placement(atom_0, atom_0.coordinate[i].x, atom_0.coordinate[i].y, atom_0.coordinate[i].z);




        double dx, dy, dz;

        dx = atom_0.coordinate[i].x - atom_0.coordinate_init[i].x;
        dy = atom_0.coordinate[i].y - atom_0.coordinate_init[i].y;
        dz = atom_0.coordinate[i].z - atom_0.coordinate_init[i].z;

        mirror_constrain_copy(atom_0.box, dx, dy, dz);// prevent the atom pass through

        double dr_square = dx * dx + dy * dy + dz * dz;
        if (dr_square > atom_0.max_nl_dr_square)
        atom_0.max_nl_dr_square = dr_square;   

    }
}








void algorithm::velvet_LJ_nl(double &dt, atom &atom_0)
{

    velvet_for_update_v_first_step(atom_0, dt);//update velocity


    if (atom_0.max_nl_dr_square * 4 > coefficient::nl_delta_square)
    {
        atom_0.max_nl_dr_square = 0;
        delete[] atom_0.neighbour_list;
        atom_0.initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = atom_0.coordinate[i].x * coefficient::cutoff_box_inv_x;
        coodinate[1] = atom_0.coordinate[i].y * coefficient::cutoff_box_inv_y;
        coodinate[2] = atom_0.coordinate[i].z * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == atom_0.boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == atom_0.boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == atom_0.boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        #pragma omp critical
        {
            atom_0.insert_in_boxes(coodinate, i);
        }
        }
        
        atom_0.initialise_neighbour_list();
        atom_0.update_neighbour_list();
        delete[] atom_0.boxes;  //free the memory
        delete[] atom_0.coordinate_init;
        atom_0.coordinate_init = new three_dim_vector[atom_0.total_num];
        for (int i = 0; i < atom_0.total_num; ++i)
        {
            atom_0.coordinate_init[i] = atom_0.coordinate[i];
        }
    }



    potential_energy::LJ_potential_nl(atom_0); //update force


    velvet_for_update_v(atom_0, dt);//update velocity

}




    void algorithm::velvet_LJ(double &dt, atom &atom_0)
    {
    velvet_for_update_v(atom_0, dt);//update velocity

    velvet_for_update_place(dt, atom_0); //update the place

    potential_energy::LJ_potential(atom_0); //update force

    velvet_for_update_v(atom_0, dt);//update velocity
    }



    void algorithm::velvet_EAM(double& dt, atom& atom_0)
    {
        velvet_for_update_v(atom_0, dt);//update velocity

        velvet_for_update_place(dt, atom_0); //update the place

        potential_energy::EAM_potential_cubicspline(atom_0); //update force

        velvet_for_update_v(atom_0, dt);//update velocity
    }




void algorithm::velvet_EAM_nl(double & dt, atom &atom_0)
{
    velvet_for_update_v_first_step(atom_0, dt);//update velocity and placement



    if (atom_0.max_nl_dr_square * 4 > coefficient::nl_delta_square)
    {
        atom_0.max_nl_dr_square = 0;
        delete[] atom_0.neighbour_list;
        atom_0.initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = atom_0.coordinate[i].x * coefficient::cutoff_box_inv_x;
        coodinate[1] = atom_0.coordinate[i].y * coefficient::cutoff_box_inv_y;
        coodinate[2] = atom_0.coordinate[i].z * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == atom_0.boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == atom_0.boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == atom_0.boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        #pragma omp critical
        {
            atom_0.insert_in_boxes(coodinate, i);
        }
        }
        
        atom_0.initialise_neighbour_list();
        atom_0.update_neighbour_list();
        delete[] atom_0.boxes;  //free the memory
        delete[] atom_0.coordinate_init;
        atom_0.coordinate_init = new three_dim_vector[atom_0.total_num];
        for (int i = 0; i < atom_0.total_num; ++i)
        {
            atom_0.coordinate_init[i] = atom_0.coordinate[i];
        }
    }

    potential_energy::EAM_nl_linear_interpolation(atom_0); //update force


    velvet_for_update_v(atom_0, dt);//update velocity
}



void algorithm::velvet_EAM_alloy(double& dt, atom& atom_0)
{
    velvet_for_update_v(atom_0, dt);//update velocity

    velvet_for_update_place(dt, atom_0); //update the place

    potential_energy::EAM_alloy_potential_cubicspline(atom_0); //update force

    velvet_for_update_v(atom_0, dt);//update velocity
}



void algorithm::velvet_EAM_alloy_nl(double & dt, atom &atom_0)
{
    velvet_for_update_v_first_step(atom_0, dt);//update velocity and placement



    if (atom_0.max_nl_dr_square * 4 > coefficient::nl_delta_square)
    {
        atom_0.max_nl_dr_square = 0;
        delete[] atom_0.neighbour_list;
        atom_0.initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = atom_0.coordinate[i].x * coefficient::cutoff_box_inv_x;
        coodinate[1] = atom_0.coordinate[i].y * coefficient::cutoff_box_inv_y;
        coodinate[2] = atom_0.coordinate[i].z * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == atom_0.boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == atom_0.boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == atom_0.boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        #pragma omp critical
        {
            atom_0.insert_in_boxes(coodinate, i);
        }
        }
        
        atom_0.initialise_neighbour_list();
        atom_0.update_neighbour_list();
        
        delete[] atom_0.boxes;  //free the memory
        delete[] atom_0.coordinate_init;
        atom_0.coordinate_init = new three_dim_vector[atom_0.total_num];
        for (int i = 0; i < atom_0.total_num; ++i)
        {
            atom_0.coordinate_init[i] = atom_0.coordinate[i];
        }
    }

    potential_energy::EAM_alloy_nl_cubic_spline(atom_0); //update force


    velvet_for_update_v(atom_0, dt);//update velocity   
}