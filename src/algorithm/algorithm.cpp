#include "../../include/algorithm/algorithm.h"
#include <iostream>
#include <algorithm>



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
        atom_0.vx[i] += 0.5 * (atom_0.fx[i] * atom_0.mass_inv[i]) * t; 
        atom_0.vy[i] += 0.5 * (atom_0.fy[i] * atom_0.mass_inv[i]) * t; 
        atom_0.vz[i] += 0.5 * (atom_0.fz[i] * atom_0.mass_inv[i]) * t; 
    }
}




//update the placement
inline void algorithm::velvet_for_update_place(double & dt, atom& atom_0)
{
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        double displacement_old[3];
        displacement_old[0] = atom_0.vx[i] * dt;
        displacement_old[1] = atom_0.vy[i] * dt;
        displacement_old[2] = atom_0.vz[i] * dt;       

        double displacement_new[3];
        compute_1dv_times_3dm(displacement_old, atom_0.basic_vectors_inv, displacement_new);

        atom_0.x[i] += displacement_new[0];
        atom_0.y[i] += displacement_new[1];
        atom_0.z[i] += displacement_new[2];
        mirror_constrain_for_placement(atom_0, atom_0.x[i], atom_0.y[i], atom_0.z[i]);
    }
}



//update the placement with neighbour list
inline void algorithm::velvet_for_update_place_with_nl(double & dt, atom& atom_0)
{
    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        double displacement_old[3];
        displacement_old[0] = atom_0.vx[i] * dt;
        displacement_old[1] = atom_0.vy[i] * dt;
        displacement_old[2] = atom_0.vz[i] * dt;  
             
        double displacement_new[3];
        compute_1dv_times_3dm(displacement_old, atom_0.basic_vectors_inv, displacement_new);


        atom_0.x[i] += displacement_new[0];
        atom_0.y[i] += displacement_new[1];
        atom_0.z[i] += displacement_new[2];
        mirror_constrain_for_placement(atom_0, atom_0.x[i], atom_0.y[i], atom_0.z[i]);

        double dr_square = (atom_0.x[i] - atom_0.x0[i]) * (atom_0.x[i] - atom_0.x0[i])
        + (atom_0.y[i] - atom_0.y0[i]) * (atom_0.y[i] - atom_0.y0[i])
        + (atom_0.z[i] - atom_0.z0[i]) * (atom_0.z[i] - atom_0.z0[i]);
        if (dr_square > atom_0.max_nl_dr_square)
        atom_0.max_nl_dr_square = dr_square;   
    }
}





void algorithm::velvet_LJ_nl(double &dt, atom &atom_0)
{

    velvet_for_update_v(atom_0, dt);//update velocity

    velvet_for_update_place_with_nl(dt, atom_0); //update the place

    if (atom_0.max_nl_dr_square * 4 > coefficient::cut_off_nl_square)
    {
        delete[] atom_0.neighbour_list;
        atom_0.initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = atom_0.x[i] * coefficient::cutoff_box_inv_x;
        coodinate[1] = atom_0.y[i] * coefficient::cutoff_box_inv_y;
        coodinate[2] = atom_0.z[i] * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == atom_0.boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == atom_0.boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == atom_0.boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        atom_0.insert_in_boxes(coodinate, i);
        }
        
        atom_0.initialise_neighbour_list();
        atom_0.update_neighbour_list();
        delete[] atom_0.boxes;  //free the memory
        atom_0.x0 = atom_0.x;
        atom_0.y0 = atom_0.y;
        atom_0.z0 = atom_0.z;
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
    velvet_for_update_v(atom_0, dt);//update velocity

    velvet_for_update_place_with_nl(dt, atom_0); //update the place

    if (atom_0.max_nl_dr_square * 4 > coefficient::cut_off_nl_square)
    {
        delete[] atom_0.neighbour_list;
        atom_0.initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < atom_0.total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = atom_0.x[i] * coefficient::cutoff_box_inv_x;
        coodinate[1] = atom_0.y[i] * coefficient::cutoff_box_inv_y;
        coodinate[2] = atom_0.z[i] * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == atom_0.boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == atom_0.boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == atom_0.boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        atom_0.insert_in_boxes(coodinate, i);
        }
        
        atom_0.initialise_neighbour_list();
        atom_0.update_neighbour_list();
        delete[] atom_0.boxes;  //free the memory
        atom_0.x0 = atom_0.x;
        atom_0.y0 = atom_0.y;
        atom_0.z0 = atom_0.z;
    }

    potential_energy::EAM_nl_linear_interpolation(atom_0); //update force


    velvet_for_update_v(atom_0, dt);//update velocity
}