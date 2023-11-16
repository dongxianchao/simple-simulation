#include "../../include/atom/atom_potential.h"
#include <iostream>
#include <random>

//this is used for inversing the matrix
std::vector<std::vector<double>> matrix_inv(std::vector<std::vector<double>> & mat)
{
    double a = mat[0][0], b = mat[0][1], c = mat[0][2];
    double d = mat[1][0], e = mat[1][1], f = mat[1][2];
    double g = mat[2][0], h = mat[2][1], i = mat[2][2];


    double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);//   determinant

    std::vector<std::vector<double>> inv(3, std::vector<double>(3));

    inv[0][0] = (e*i - f*h) / det;
    inv[0][1] = (c*h - b*i) / det;
    inv[0][2] = (b*f - c*e) / det;
    inv[1][0] = (f*g - d*i) / det;
    inv[1][1] = (a*i - c*g) / det;
    inv[1][2] = (c*d - a*f) / det;
    inv[2][0] = (d*h - e*g) / det;
    inv[2][1] = (b*g - a*h) / det;
    inv[2][2] = (a*e - b*d) / det;

    return inv;
}

//this is used for store datas
struct vec
{
    double mold, x, y, z, force;
};


void atom::update_the_virial_tensor(double& Fx, double& Fy, double& Fz, double& rx, double& ry, double& rz, int& index)
{
    virial_tensor[index][0] += rx * Fx;
    virial_tensor[index][1] += rx * Fy;
    virial_tensor[index][2] += rx * Fz;
    virial_tensor[index][3] += ry * Fx;
    virial_tensor[index][4] += ry * Fy;
    virial_tensor[index][5] += ry * Fz;
    virial_tensor[index][6] += rz * Fx;
    virial_tensor[index][7] += rz * Fy;
    virial_tensor[index][8] += rz * Fz;
}

void atom::clean_the_virial_tensor()
{
    for (int i = 0; i < total_num; ++i)
    std::fill(virial_tensor[i].begin(), virial_tensor[i].end(), 0);
}


atom::atom(std::vector<int> numcells, int numbers,
    std::vector<double> mass_0,
    std::vector<double> x_0,
    std::vector<double> y_0,
    std::vector<double> z_0, cell& cell_0)
    {
        primitive_cell_num = numcells;
        cell_num = numbers;
        total_num = cell_num * numcells[0] * numcells[1] * numcells[2];
        coordinate = new three_dim_vector[total_num];
        velocity = new three_dim_vector[total_num];
        coordinate_init = new three_dim_vector[total_num];
        force = new three_dim_vector[total_num];
        mass.resize(total_num, 0);
        count.resize(total_num, 0);
        mass_inv.resize(total_num, 0);
        pe.resize(total_num, 0);
        type.resize(total_num, 0);
        rho.resize(total_num, 0);   // allocate memory
        virial_tensor.reserve(total_num);

        //for (int i = 0; i < total_num; ++i)
        //std::fill(virial_tensor[i], virial_tensor[i] + 9, 0);

        cell_in = &cell_0;

        if_othogonal = false;

        if (cell_0.matrix == std::vector<std::vector<double>>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}})
        {
            if_othogonal = true;
        }// set the condition



        box[0] = numcells[0] * cell_0.lattice_x;
        box[1] = numcells[1] * cell_0.lattice_y;
        box[2] = numcells[2] * cell_0.lattice_z;
        box[3] = 0.5 * box[0];
        box[4] = 0.5 * box[1];
        box[5] = 0.5 * box[2];  //initialise box




        basic_vectors = cell_0.matrix;
        basic_vectors_inv = matrix_inv(basic_vectors);//    get the matrix




        int n = 0;
        for (int ix = 0; ix < numcells[0]; ++ix)
        {
            for (int iy = 0; iy < numcells[1]; ++iy)
            {
                for (int iz = 0; iz < numcells[2]; ++iz)
                {
                    for (int i = 0; i < numbers; ++i)
                    {
                        coordinate[n].x = (ix + x_0[i]) * cell_0.lattice_x;
                        coordinate[n].y = (iy + y_0[i]) * cell_0.lattice_y;
                        coordinate[n].z = (iz + z_0[i]) * cell_0.lattice_z;
                        mass[n] = mass_0[i];
                        mass_inv[n] = 1 / mass_0[i];
                        ++n;
                    }
                }
            }
        }// used for insert every atom in their own place.
    }


    double atom::kinetic_energy() //compute the kinetic energy
    {
        double energy = 0;

        for (int i = 0; i < total_num; ++i)
        {
            energy += 0.5 * (velocity[i].x * velocity[i].x + velocity[i].y * velocity[i].y + velocity[i].z * velocity[i].z) * mass[i];
        }
        kineticenergy = energy;
        return energy;
    }

    double atom::find_temperature()
    {
        T = 2 * kineticenergy / (3 * coefficient::k_B * total_num);
        return T;
    }


    void atom::initialise_velocity(double T0)
    {
        std::vector<double> centerofmassvelocity = {0, 0, 0};
        double total_mass = 0;
        for (int i = 0; i < total_num; ++i)        // add the mass up
        {
            total_mass += mass[i];
        }


        for (int n = 0; n < total_num; ++n)  // randomize the velocity and compute the center velocity
        {
            velocity[n].x = -1 + (rand() * 2.0) / RAND_MAX;
            velocity[n].y = -1 + (rand() * 2.0) / RAND_MAX;
            velocity[n].z = -1 + (rand() * 2.0) / RAND_MAX;
            centerofmassvelocity[0] += velocity[n].x * mass[n];
            centerofmassvelocity[1] += velocity[n].y * mass[n];
            centerofmassvelocity[2] += velocity[n].z * mass[n];
        }
        centerofmassvelocity[0] /= total_mass;
        centerofmassvelocity[1] /= total_mass;
        centerofmassvelocity[2] /= total_mass;



        for (int i = 0; i < total_num; ++i) //compute the relative velocity
        {
            velocity[i].x -= centerofmassvelocity[0];
            velocity[i].y -= centerofmassvelocity[1];
            velocity[i].z -= centerofmassvelocity[2];
        }



        double energy = kinetic_energy();// calculate the temperature
        double t = 2 * energy / (3 * coefficient::k_B * total_num);

        double scalar = sqrt(T0 / t);       //find the velocity scalar



        for (int i = 0; i < total_num; ++i)     //update the velocity
        {
            velocity[i].x *= scalar;
            velocity[i].y *= scalar;
            velocity[i].z *= scalar;
        }
    }





// just as the name says
void compute_1dv_times_3dm(double* arr, std::vector<std::vector<double>> &mtx, double (&arr_out)[3]) {
    for (int i = 0; i < 3; ++i) {
        arr_out[i] = arr[0] * mtx[0][i] + arr[1] * mtx[1][i] + arr[2] * mtx[2][i];
    }
}



//do some auxilary work
inline void mirror_constrain_once(double& len,double& half_len, double &xij)
{
    if (xij > half_len)
    xij -= len;
    else if (xij < -half_len)
    xij += len;
}

void mirror_constrain(double* b, double &xij, double &yij, double &zij)
// to apply the mirror constrain
{
    mirror_constrain_once(b[0], b[3], xij);
    mirror_constrain_once(b[1], b[4], yij);
    mirror_constrain_once(b[2], b[5], zij);
}



// it's the classical LJ method which will calculate the all atoms
void potential_energy::LJ_potential(atom & atom_0)
{
    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.force[n].x = atom_0.force[n].y = atom_0.force[n].z = atom_0.pe[n] = 0;
    }

    for (int i = 0; i < atom_0.total_num - 1; ++i)  //calculate the force
    {
        for (int j = i + 1; j < atom_0.total_num; ++j)
        {
            double xij = atom_0.coordinate[j].x - atom_0.coordinate[i].x;
            double yij = atom_0.coordinate[j].y - atom_0.coordinate[i].y;
            double zij = atom_0.coordinate[j].z - atom_0.coordinate[i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            double old_r[3] = {xij, yij, zij};
            double new_r[3];
            compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
            xij = new_r[0];
            yij = new_r[1];
            zij = new_r[2];
            double r2 = xij*xij + yij*yij + zij*zij;


            if (r2 > coefficient::cutoffSquare)
            continue;               //set the cutoff


            double r2inv = 1.0 / r2;
            const double r4inv = r2inv * r2inv;
            const double r6inv = r2inv * r4inv;
            const double r8inv = r4inv * r4inv;
            const double r12inv = r4inv * r8inv;
            const double r14inv = r6inv * r8inv;
            const double f_ij = coefficient::e24s6 * r8inv - coefficient::e48s12 * r14inv;
            atom_0.pe[i] += coefficient::e4s12 * r12inv - coefficient::e4s6 * r6inv;

            double Fx, Fy, Fz;
            Fx = f_ij * xij;
            Fy = f_ij * yij;
            Fz = f_ij * zij;
            atom_0.force[i].x += Fx;
            atom_0.force[j].x -= Fx;
            atom_0.force[i].y += Fy;
            atom_0.force[j].y -= Fy;
            atom_0.force[i].z += Fz;
            atom_0.force[j].z -= Fz;

            if (atom_0.virial_switch)
            {
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, xij, yij, zij, i);
            }

        }
    }
}

// calculate the energy
double atom::potential_energy()
{
    double energy = 0;

    for (int i = 0; i < total_num; ++i)
    energy += pe[i];
    potentialenergy = energy;
    return energy;
}
//this either
double atom::total_energy()
{
    totalenergy = kineticenergy + potentialenergy;
    return totalenergy;
}

void atom::calculate_virial_pressure()
{
    //initialise
    for (int i = 0; i < 9; ++i)
    virial_press_tensor_tot[i] = 0;

    for (int i = 0; i < total_num; ++i)
    {
        for(int j = 0; j < 9; ++j)
        virial_press_tensor_tot[j] -= virial_tensor[i][j];
    }

    double vxx, volume; //this is the second part of the equation

    double a = box[0] * basic_vectors[0][0], b = box[0] * basic_vectors[0][1], c = box[0] * basic_vectors[0][2];
    double d = box[1] * basic_vectors[1][0], e = box[1] * basic_vectors[1][1], f = box[1] * basic_vectors[1][2];
    double g = box[2] * basic_vectors[2][0], h = box[2] * basic_vectors[2][1], i = box[2] * basic_vectors[2][2];

    volume = abs(a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g));//   determinant

    double volume_inv = 1 / volume;

    vxx = total_num * coefficient::k_B * T;

    virial_press_tensor_tot[0] += vxx;
    virial_press_tensor_tot[3] += vxx;
    virial_press_tensor_tot[6] += vxx;

    for (int i = 0; i < 9; ++i)
    virial_press_tensor_tot[i] *= volume_inv;

}

//input miu matrix and vector and output the outcome
void util_calculator_for_rescale(double* input_miu, three_dim_vector& v, double* output)
{
    output[0] = input_miu[0] * v.x + input_miu[1] * v.y + input_miu[2] * v.z;
    output[1] = input_miu[3] * v.x + input_miu[4] * v.y + input_miu[5] * v.z;
    output[2] = input_miu[6] * v.x + input_miu[7] * v.y + input_miu[8] * v.z;
}

void atom::rescale_the_position(double* miu)
{
    max_nl_dr_square = 0;
    for (int i = 0; i < total_num; ++i)
    {
        double final_coordinate[3];
        double raw[3] = {coordinate[i].x, coordinate[i].y, coordinate[i].z};//this coordinate hasn't been transformed
        double init_coordinate[3];

        compute_1dv_times_3dm(raw, basic_vectors, init_coordinate);
        util_calculator_for_rescale(miu, coordinate[i], final_coordinate);

        double dx, dy, dz;
        dx = final_coordinate[0] - init_coordinate[0];
        dy = final_coordinate[1] - init_coordinate[1];
        dz = final_coordinate[2] - init_coordinate[2];

        double final_coordinate_transform[3];

        compute_1dv_times_3dm(final_coordinate, basic_vectors_inv, final_coordinate_transform);

        coordinate[i].x = final_coordinate_transform[0];
        coordinate[i].y = final_coordinate_transform[1];
        coordinate[i].z = final_coordinate_transform[2];

        double dr_square = dx * dx + dy * dy + dz * dz;
        if (dr_square > max_nl_dr_square)
        max_nl_dr_square = dr_square;
    }// this part is used for update position and decide if need to update the nl

    cell_in->matrix[0][0] *= cell_in->lattice_x;
    cell_in->matrix[0][1] *= cell_in->lattice_x;
    cell_in->matrix[0][2] *= cell_in->lattice_x;
    cell_in->matrix[1][0] *= cell_in->lattice_y;
    cell_in->matrix[1][1] *= cell_in->lattice_y;
    cell_in->matrix[1][2] *= cell_in->lattice_y;
    cell_in->matrix[2][0] *= cell_in->lattice_z;
    cell_in->matrix[2][1] *= cell_in->lattice_z;
    cell_in->matrix[2][2] *= cell_in->lattice_z;//recalculate the box, lattice and the vector

    std::vector<std::vector<double>> new_mtx(3, std::vector<double>(3));
    new_mtx[0][0] = miu[0] * cell_in->matrix[0][0] + miu[1] * cell_in->matrix[0][1] + miu[2] * cell_in->matrix[0][2];
    new_mtx[1][0] = miu[0] * cell_in->matrix[1][0] + miu[1] * cell_in->matrix[1][1] + miu[2] * cell_in->matrix[1][2];
    new_mtx[2][0] = miu[0] * cell_in->matrix[2][0] + miu[1] * cell_in->matrix[2][1] + miu[2] * cell_in->matrix[2][2];
    new_mtx[0][1] = miu[3] * cell_in->matrix[0][0] + miu[4] * cell_in->matrix[0][1] + miu[5] * cell_in->matrix[0][2];
    new_mtx[1][1] = miu[3] * cell_in->matrix[1][0] + miu[4] * cell_in->matrix[1][1] + miu[5] * cell_in->matrix[1][2];
    new_mtx[2][1] = miu[3] * cell_in->matrix[2][0] + miu[4] * cell_in->matrix[2][1] + miu[5] * cell_in->matrix[2][2];
    new_mtx[0][2] = miu[6] * cell_in->matrix[0][0] + miu[7] * cell_in->matrix[0][1] + miu[8] * cell_in->matrix[0][2];
    new_mtx[1][2] = miu[6] * cell_in->matrix[1][0] + miu[7] * cell_in->matrix[1][1] + miu[8] * cell_in->matrix[1][2];
    new_mtx[2][2] = miu[6] * cell_in->matrix[2][0] + miu[7] * cell_in->matrix[2][1] + miu[8] * cell_in->matrix[2][2];
    //it's right because of the transpose


    cell_in->lattice_x = sqrt(new_mtx[0][0] * new_mtx[0][0] + new_mtx[0][1] * new_mtx[0][1] + new_mtx[0][2] * new_mtx[0][2]);
    cell_in->lattice_y = sqrt(new_mtx[1][0] * new_mtx[1][0] + new_mtx[1][1] * new_mtx[1][1] + new_mtx[1][2] * new_mtx[1][2]);
    cell_in->lattice_z = sqrt(new_mtx[2][0] * new_mtx[2][0] + new_mtx[2][1] * new_mtx[2][1] + new_mtx[2][2] * new_mtx[2][2]);


    box[0] = cell_in->lattice_x * primitive_cell_num[0];
    box[1] = cell_in->lattice_y * primitive_cell_num[1];
    box[2] = cell_in->lattice_z * primitive_cell_num[2];
    box[3] = 0.5 * box[0];
    box[4] = 0.5 * box[1];
    box[5] = 0.5 * box[2];

    for (auto& v : new_mtx)
    {
        double magnitude = 0.0;
        for (double val : v) {
            magnitude += val * val;
        }
    magnitude = std::sqrt(magnitude);

    v[0] /= magnitude;
    v[1] /= magnitude;
    v[2] /= magnitude;
    }//normalize the matrix

    cell_in->matrix = new_mtx;
    basic_vectors = new_mtx;
    basic_vectors_inv = matrix_inv(new_mtx);

    initialise_boxes_xyz_maxnums(*cell_in);
    
    if (max_nl_dr_square * 4 > coefficient::nl_delta_square)
    {
        max_nl_dr_square = 0;
        delete[] neighbour_list;
        initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = coordinate[i].x * coefficient::cutoff_box_inv_x;
        coodinate[1] = coordinate[i].y * coefficient::cutoff_box_inv_y;
        coodinate[2] = coordinate[i].z * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        #pragma omp critical
        {
            insert_in_boxes(coodinate, i);
        }
        }
        
        initialise_neighbour_list();
        update_neighbour_list();
        delete[] boxes;  //free the memory
        delete[] coordinate_init;
        coordinate_init = new three_dim_vector[total_num];
        for (int i = 0; i < total_num; ++i)
        {
            coordinate_init[i] = coordinate[i];
        }
    }
}


void atom::rescale_the_position(double miuxx, double miuyy, double miuzz)
{
    max_nl_dr_square = 0;
    for (int i = 0; i < total_num; ++i)
    {
        double dx, dy, dz;
        dx = (miuxx - 1) * coordinate[i].x;
        dy = (miuyy - 1) * coordinate[i].y;
        dz = (miuzz - 1) * coordinate[i].z;
        coordinate[i].x *= miuxx;
        coordinate[i].y *= miuyy;
        coordinate[i].z *= miuzz;

        double dr_square = dx * dx + dy * dy + dz * dz;
        if (dr_square > max_nl_dr_square)
        max_nl_dr_square = dr_square;
    }// this part is used for update position and decide if need to update the nl

    box[0] *= miuxx;
    box[1] *= miuyy;
    box[2] *= miuzz;
    box[3] = 0.5 * box[0];
    box[4] = 0.5 * box[1];
    box[5] = 0.5 * box[2];
    cell_in->lattice_x *= miuxx;
    cell_in->lattice_y *= miuyy;
    cell_in->lattice_z *= miuzz;


    initialise_boxes_xyz_maxnums(*cell_in);
    
    if (max_nl_dr_square * 4 > coefficient::nl_delta_square)
    {
        max_nl_dr_square = 0;
        delete[] neighbour_list;
        initialise_boxes();

        #pragma omp parallel for
        for (int i = 0; i < total_num; ++i)
        {
        int coodinate[3];// used for insert in box
        coodinate[0] = coordinate[i].x * coefficient::cutoff_box_inv_x;
        coodinate[1] = coordinate[i].y * coefficient::cutoff_box_inv_y;
        coodinate[2] = coordinate[i].z * coefficient::cutoff_box_inv_z;

        if (coodinate[0] == boxes_x) coodinate[0] -= 1;
        if (coodinate[1] == boxes_y) coodinate[1] -= 1;
        if (coodinate[2] == boxes_z) coodinate[2] -= 1;//    apply the method of mirror

        #pragma omp critical
        {
            insert_in_boxes(coodinate, i);
        }
        }
        
        initialise_neighbour_list();
        update_neighbour_list();
        delete[] boxes;  //free the memory
        delete[] coordinate_init;
        coordinate_init = new three_dim_vector[total_num];
        for (int i = 0; i < total_num; ++i)
        {
            coordinate_init[i] = coordinate[i];
        }
    }
}

void atom::update_energy()
{
    kinetic_energy();
    find_temperature();
    potential_energy();
    total_energy();
}



inline double vector_dot(std::vector<double> & v1, std::vector<double> & v2)//  as the name said
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline std::vector<double> vector_cross(std::vector<double> & v1, std::vector<double> & v2)//as the name said
{
    std::vector<double> v3 = {  v1[1] * v2[2] - v1[2] * v2[1],
                                v1[2] * v2[0] - v1[0] * v2[2],
                                v1[0] * v2[1] - v1[1] * v2[0]};
    return v3;
}

inline std::vector<double> opposite_vector(std::vector<double> & v1)// used for an opposite vector
{
    std::vector<double> _v = {-v1[0], -v1[1], -v1[2]};
    return _v;
}

double mold_vector(std::vector<double> v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


void atom::initialise_boxes_xyz_maxnums(cell& cell_0)
{
    double dotxy, dotyz, dotxz;
    dotxy = vector_dot(cell_0.matrix[0], cell_0.matrix[1]);
    dotyz = vector_dot(cell_0.matrix[1], cell_0.matrix[2]);
    dotxz = vector_dot(cell_0.matrix[0], cell_0.matrix[2]);
    if (dotxy != 0)
    {
        double cos_xoy_minus_half_pi = mold_vector(vector_cross(cell_0.matrix[0], cell_0.matrix[1]));
        coefficient::cutoff_box_x = ((coefficient::cutoff + 1) / cos_xoy_minus_half_pi > coefficient::cutoff_box_x)? (coefficient::cutoff + 1) / cos_xoy_minus_half_pi : coefficient::cutoff_box_x;
        coefficient::cutoff_box_y = ((coefficient::cutoff + 1) / cos_xoy_minus_half_pi > coefficient::cutoff_box_y)? (coefficient::cutoff + 1) / cos_xoy_minus_half_pi : coefficient::cutoff_box_y;
    }
   if (dotxz != 0)
    {
        double cos_xoz_minus_half_pi = mold_vector(vector_cross(cell_0.matrix[0], cell_0.matrix[2]));
        coefficient::cutoff_box_x = ((coefficient::cutoff + 1) / cos_xoz_minus_half_pi > coefficient::cutoff_box_x)? (coefficient::cutoff + 1) / cos_xoz_minus_half_pi : coefficient::cutoff_box_x;
        coefficient::cutoff_box_z = ((coefficient::cutoff + 1) / cos_xoz_minus_half_pi > coefficient::cutoff_box_z)? (coefficient::cutoff + 1) / cos_xoz_minus_half_pi : coefficient::cutoff_box_z;
    }
   if (dotyz != 0)
    {
        double cos_yoz_minus_half_pi = mold_vector(vector_cross(cell_0.matrix[1], cell_0.matrix[2]));
        coefficient::cutoff_box_y = ((coefficient::cutoff + 1) / cos_yoz_minus_half_pi > coefficient::cutoff_box_y)? (coefficient::cutoff + 1) / cos_yoz_minus_half_pi : coefficient::cutoff_box_y;
        coefficient::cutoff_box_z = ((coefficient::cutoff + 1) / cos_yoz_minus_half_pi > coefficient::cutoff_box_z)? (coefficient::cutoff + 1) / cos_yoz_minus_half_pi : coefficient::cutoff_box_z;
    }


    boxes_x = box[0] / coefficient::cutoff_box_x;
    boxes_y = box[1] / coefficient::cutoff_box_y;
    boxes_z = box[2] / coefficient::cutoff_box_z;
    coefficient::cutoff_box_inv_x = 1 / coefficient::cutoff_box_x;
    coefficient::cutoff_box_inv_y = 1 / coefficient::cutoff_box_y;
    coefficient::cutoff_box_inv_z = 1 / coefficient::cutoff_box_z;

}



void atom::initialise_boxes()
{
    boxes = new int[boxes_x * boxes_y * boxes_z * coefficient::max_num_in_boxes]();// create the ptr
}








int* atom::get_place(double& x, double& y, double& z, atom& atom_0)
//  to get the boxes_index of the specific coodinate
{
    int* out = new int[3];
    out[0] = x * coefficient::cutoff_box_inv_x;
    out[1] = y * coefficient::cutoff_box_inv_y;
    out[2] = z * coefficient::cutoff_box_inv_z;
    if (out[0] == atom_0.boxes_x)
    out[0] -= 1;
    if (out[1] == atom_0.boxes_y)
    out[1] -= 1;
    if (out[2] == atom_0.boxes_z)
    out[2] -= 1;        //don't delete this if code
    return out;
}

void atom::insert_in_boxes(int*& place, int index)
{
    int prefix = (place[0] * boxes_y * boxes_z + place[1] * boxes_z + place[2])
    * coefficient::max_num_in_boxes;
    if (boxes[prefix] >= coefficient::max_num_in_boxes - 1)
    throw std::out_of_range("numbers in box out of range");// for debuging
    boxes[prefix] += 1;
    boxes[boxes[prefix] + prefix] = index;
}


void atom::insert_in_boxes(int (&place)[3], int index)
{
    int prefix = (place[0] * boxes_y * boxes_z + place[1] * boxes_z + place[2])
    * coefficient::max_num_in_boxes;
    if (boxes[prefix] >= coefficient::max_num_in_boxes - 1)
    throw std::out_of_range("numbers in box out of range");// for debugging
    boxes[prefix] += 1;
    boxes[boxes[prefix] + prefix] = index;
}// used for insert in boxes



void check_if_on_the_edge(int& x, int& y, int& z, atom& atom_0)
{
    if (x == -1 || x == atom_0.boxes_x)
    x = (x == -1)? atom_0.boxes_x -1 : 0;
    if (y == -1 || y == atom_0.boxes_y)
    y = (y == -1)? atom_0.boxes_y -1 : 0;
    if (z == -1 || z == atom_0.boxes_z)
    z = (z == -1)? atom_0.boxes_z -1 : 0;
}



void potential_energy::LJ_potential_nl(atom & atom_0)
{

    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();

    delete[] atom_0.force;
    atom_0.force = new three_dim_vector[atom_0.total_num]();
    //#pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.pe[n] = 0;
    }// initialise the data

    #pragma omp parallel for
    for (int index_i = 0; index_i < atom_0.total_num; ++index_i)  //calculate the force
    {
        int prefix = index_i * coefficient::max_num_in_nl;
        double local_force_x = 0, local_force_y = 0, local_force_z = 0;
        for(int s = 1; s <= atom_0.neighbour_list[prefix]; ++s)
        {
            int index_j = atom_0.neighbour_list[prefix + s];
            double xij = atom_0.coordinate[index_j].x - atom_0.coordinate[index_i].x;
            double yij = atom_0.coordinate[index_j].y - atom_0.coordinate[index_i].y;
            double zij = atom_0.coordinate[index_j].z - atom_0.coordinate[index_i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            double old_r[3] = {xij, yij, zij};
            double new_r[3];
            compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
            xij = new_r[0];
            yij = new_r[1];
            zij = new_r[2];
            double r2 = xij*xij + yij*yij + zij*zij;

            if (r2 > coefficient::cutoffSquare)
            continue;


            double r2inv = 1.0 / r2;
            const double r4inv = r2inv * r2inv;
            const double r6inv = r2inv * r4inv;
            const double r8inv = r4inv * r4inv;
            const double r12inv = r4inv * r8inv;
            const double r14inv = r6inv * r8inv;
            const double f_ij = coefficient::e24s6 * r8inv - coefficient::e48s12 * r14inv;
            #pragma omp atomic
            atom_0.pe[index_i] += 0.5 * (coefficient::e4s12 * r12inv - coefficient::e4s6 * r6inv);
            #pragma omp atomic
            atom_0.pe[index_i] += 0.5 * (coefficient::e4s12 * r12inv - coefficient::e4s6 * r6inv);


            double Fx, Fy, Fz;
            Fx = f_ij * xij;
            Fz = f_ij * yij;
            Fz = f_ij * zij;

            local_force_x += Fx;
            local_force_y += Fy;
            local_force_z += Fz;

            #pragma omp atomic
            atom_0.force[index_j].x -= f_ij * xij;
            #pragma omp atomic
            atom_0.force[index_j].y -= f_ij * yij;
            #pragma omp atomic
            atom_0.force[index_j].z -= f_ij * zij;

            if (atom_0.virial_switch)
            {
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, xij, yij, zij, index_i);
            }


        }
        #pragma omp atomic
        atom_0.force[index_i].x += local_force_x;
        #pragma omp atomic
        atom_0.force[index_i].y += local_force_y;
        #pragma omp atomic
        atom_0.force[index_i].z += local_force_z;
    }
}



void atom::initialise_neighbour_list()// used to initialise nl
{
    neighbour_list = new int[coefficient::max_num_in_nl * total_num]();
}

void atom::update_neighbour_list()
{
    #pragma omp parallel for
    for (int x = 0; x < boxes_x; ++x)
    {
        for (int y = 0; y < boxes_y; ++y)
        {
            for (int z = 0; z < boxes_z; ++z)
            {
                int prefix = (x * boxes_y * boxes_z + y * boxes_z + z ) * coefficient::max_num_in_boxes;
                int num = boxes[prefix];
                for (int i = 1; i <= num - 1; ++i)
                {
                    int index_i = boxes[prefix + i];
                    for(int j = i + 1; j <= num; ++j)
                    {
                        int index_j = boxes[prefix + j];
                        double xij = coordinate[index_j].x - coordinate[index_i].x;
                        double yij = coordinate[index_j].y - coordinate[index_i].y;
                        double zij = coordinate[index_j].z - coordinate[index_i].z;
                        if (!if_othogonal)
                        {   
                            double old_r[3] = {xij, yij, zij};
                            double new_r[3];
                            compute_1dv_times_3dm(old_r, basic_vectors, new_r);
                            xij = new_r[0];
                            yij = new_r[1];
                            zij = new_r[2];
                        }
                        double r2 = xij*xij + yij*yij + zij*zij;
                        if (r2 > coefficient::cut_off_nl_square)
                            continue;// set the cutoff

                        int pre = coefficient::max_num_in_nl * index_i;
                        if (neighbour_list[pre] >= coefficient::max_num_in_nl - 1)
                        throw std::out_of_range("neighbour list out of range");
                        neighbour_list[pre] += 1;
                        neighbour_list[neighbour_list[pre] + pre] = index_j;
                    }
                }
            }
        }
    }// the first cell
  
    #pragma omp parallel for
    for (int index_i = 0; index_i < total_num; ++index_i)  //   for every atom
    {
        int place[3];
        place[0] = coordinate[index_i].x * coefficient::cutoff_box_inv_x;
        place[1] = coordinate[index_i].y * coefficient::cutoff_box_inv_y;
        place[2] = coordinate[index_i].z * coefficient::cutoff_box_inv_z;

        if (place[0] == boxes_x) place[0] -= 1;
        if (place[1] == boxes_y) place[1] -= 1;
        if (place[2] == boxes_z) place[2] -= 1;//    the same function as the get_place function

        int prefix = (place[0] * boxes_y * boxes_z + place[1] * boxes_z + place[2]) * coefficient::max_num_in_boxes;
        int num = boxes[prefix];
        int x_index = place[0] + 1;
        int y_index = place[1];
        int z_index = place[2];
        if (x_index == -1 || x_index == boxes_x)
        x_index = (x_index == -1)? boxes_x -1 : 0;
        if (y_index == -1 || y_index == boxes_y)
        y_index = (y_index == -1)? boxes_y -1 : 0;
        if (z_index == -1 || z_index == boxes_z)
        z_index = (z_index == -1)? boxes_z -1 : 0;

        prefix = (x_index * boxes_y * boxes_z + y_index * boxes_z + z_index) * coefficient::max_num_in_boxes;
        num = boxes[prefix];
        for (int a = 1; a <= num; ++a)
        {
            int index_j = boxes[prefix + a];// dont need to prevent the same atoms
            double xij = coordinate[index_j].x - coordinate[index_i].x;
            double yij = coordinate[index_j].y - coordinate[index_i].y;
            double zij = coordinate[index_j].z - coordinate[index_i].z;
            mirror_constrain(box, xij, yij, zij);
            if (!if_othogonal)
            {
                double old_r[3] = {xij, yij, zij};
                double new_r[3];
                compute_1dv_times_3dm(old_r, basic_vectors, new_r);
                xij = new_r[0];
                yij = new_r[1];
                zij = new_r[2];
            }
            double r2 = xij*xij + yij*yij + zij*zij;
            if (r2 > coefficient::cut_off_nl_square)
                continue;// set the cutoff

            int pre = coefficient::max_num_in_nl * index_i;
            if (neighbour_list[pre] >= coefficient::max_num_in_nl - 1)
            throw std::out_of_range("neighbour list out of range");
            neighbour_list[pre] += 1;
            neighbour_list[neighbour_list[pre] + pre] = index_j;
        }// the second cell



        for (int x = -1; x < 2; ++x)
        {
            x_index = place[0] + x;
            y_index = place[1] + 1;
            z_index = place[2];
            if (x_index == -1 || x_index == boxes_x)
            x_index = (x_index == -1)? boxes_x -1 : 0;
            if (y_index == -1 || y_index == boxes_y)
            y_index = (y_index == -1)? boxes_y -1 : 0;
            if (z_index == -1 || z_index == boxes_z)
            z_index = (z_index == -1)? boxes_z -1 : 0;

            prefix = (x_index * boxes_y * boxes_z + y_index * boxes_z + z_index) * coefficient::max_num_in_boxes;
            num = boxes[prefix];
            for (int a = 1; a <= num; ++a)
            {
                int index_j = boxes[prefix + a];
                double xij = coordinate[index_j].x - coordinate[index_i].x;
                double yij = coordinate[index_j].y - coordinate[index_i].y;
                double zij = coordinate[index_j].z - coordinate[index_i].z;
                mirror_constrain(box, xij, yij, zij);
                if (!if_othogonal)
                {
                    double old_r[3] = {xij, yij, zij};
                    double new_r[3];
                    compute_1dv_times_3dm(old_r, basic_vectors, new_r);
                    xij = new_r[0];
                    yij = new_r[1];
                    zij = new_r[2];
                }
                double r2 = xij*xij + yij*yij + zij*zij;
                if (r2 > coefficient::cut_off_nl_square)
                    continue;// set the cutoff

                int pre = coefficient::max_num_in_nl * index_i;
                if (neighbour_list[pre] >= coefficient::max_num_in_nl - 1)
                throw std::out_of_range("neighbour list out of range");
                neighbour_list[pre] += 1;
                neighbour_list[neighbour_list[pre] + pre] = index_j;
            }

        }

        z_index += 1;// for the top cells
        for (int x = -1; x < 2; ++x)
        {
            for(int y = -1; y < 2; ++y)
            {
                x_index = place[0] + x;
                y_index = place[1] + y;
                if (x_index == -1 || x_index == boxes_x)
                x_index = (x_index == -1)? boxes_x -1 : 0;
                if (y_index == -1 || y_index == boxes_y)
                y_index = (y_index == -1)? boxes_y -1 : 0;
                if (z_index == -1 || z_index == boxes_z)
                z_index = (z_index == -1)? boxes_z -1 : 0;

                prefix = (x_index * boxes_y * boxes_z + y_index * boxes_z + z_index) * coefficient::max_num_in_boxes;
                num = boxes[prefix];
                for (int a = 1; a <= num; ++a)
                {
                    int index_j = boxes[prefix + a];
                    double xij = coordinate[index_j].x - coordinate[index_i].x;
                    double yij = coordinate[index_j].y - coordinate[index_i].y;
                    double zij = coordinate[index_j].z - coordinate[index_i].z;
                    mirror_constrain(box, xij, yij, zij);
                    if (!if_othogonal)
                    {
                        double old_r[3] = {xij, yij, zij};
                        double new_r[3];
                        compute_1dv_times_3dm(old_r, basic_vectors, new_r);
                        xij = new_r[0];
                        yij = new_r[1];
                        zij = new_r[2];
                    }
                    double r2 = xij*xij + yij*yij + zij*zij;
                    if (r2 > coefficient::cut_off_nl_square)
                        continue;// set the cutoff

                    int pre = coefficient::max_num_in_nl * index_i;
                    if (neighbour_list[pre] >= coefficient::max_num_in_nl - 1)
                    throw std::out_of_range("neighbour list out of range");
                    neighbour_list[pre] += 1;
                    neighbour_list[neighbour_list[pre] + pre] = index_j;
                }
            }
        }
    }
    for (int i = 0; i < total_num; ++i)
    {
        int prefix = coefficient::max_num_in_nl * i;
        int num = neighbour_list[prefix];
        if (num + 1 > coefficient::max_nums_in_cutoff)
        {
            coefficient::max_nums_in_cutoff = num + 1;
        }
    }
}





double potential_energy::EAM_drho = 0.0;
double potential_energy::EAM_dr = 0.0;
double potential_energy::EAM_drho_inv = 0.0;
double potential_energy::EAM_dr_inv = 0.0;
int potential_energy::EAM_Nrho = 0;
int potential_energy::EAM_Nr = 0;
double* potential_energy::EAM_F_rho;
double* potential_energy::EAM_Z_r;
double* potential_energy::EAM_rho_r;
double* potential_energy::EAM_cubic_spline_coefficient_r_Zr;
double* potential_energy::EAM_cubic_spline_coefficient_rho_Frho;
double* potential_energy::EAM_cubic_spline_coefficient_r_rhor;
double** potential_energy::EAM_alloy_F_rho;
double** potential_energy::EAM_alloy_rfhi_r;
double** potential_energy::EAM_alloy_rho_r;
double** potential_energy::EAM_alloy_cubic_spline_coefficient_rfhi_r;
double** potential_energy::EAM_alloy_cubic_spline_coefficient_rho_r;
double** potential_energy::EAM_alloy_cubic_spline_coefficient_F_rho;// the parameters of cubic spline

void potential_energy::get_the_EAM_data(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    file >> EAM_Nrho >> EAM_drho >> EAM_Nr >> EAM_dr >> coefficient::cutoff;
    coefficient::cutoff = (coefficient::cutoff < EAM_dr * (EAM_Nr - 1)) ? coefficient::cutoff : EAM_dr * (EAM_Nr - 1);


    EAM_F_rho = new double[EAM_Nrho];
    EAM_Z_r = new double[EAM_Nr];
    EAM_rho_r = new double[EAM_Nr];

    for (int i = 0; i < EAM_Nrho; ++i)
    file >> EAM_F_rho[i];

    for (int i = 0; i < EAM_Nr; ++i)
    file >> EAM_Z_r[i];

    for (int i = 0; i < EAM_Nr; ++i)
    file >> EAM_rho_r[i];


    file.close();

    EAM_dr_inv = 1 / EAM_dr;
    EAM_drho_inv = 1 / EAM_drho;
    coefficient::cutoff_box_x = coefficient::cutoff_box_y = coefficient::cutoff_box_z = coefficient::cutoff + 1;
    coefficient::cut_off_nl = coefficient::cutoff + coefficient::nl_delta;
    coefficient::cut_off_nl_square = coefficient::cut_off_nl * coefficient::cut_off_nl;
    coefficient::cutoffSquare = coefficient::cutoff * coefficient::cutoff;
}


int potential_energy::atom_types_nums = 0;
int potential_energy::total_pair_nums = 0;//    declare the variant
void potential_energy::get_the_EAM_alloy_data(std::ifstream& file, atom& atom_0, bool random, int type, double percentage)
{
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);//ignore the first 3 lines
    file >> atom_types_nums;//   read the type nums
    std::getline(file, line);

    file >> EAM_Nrho >> EAM_drho >> EAM_Nr >> EAM_dr >> coefficient::cutoff;//read the data
    coefficient::cutoff = (coefficient::cutoff < EAM_dr * (EAM_Nr - 1))? coefficient::cutoff : EAM_dr * (EAM_Nr - 1);


    EAM_alloy_F_rho = new double* [atom_types_nums];
    EAM_alloy_cubic_spline_coefficient_F_rho = new double* [atom_types_nums];
    EAM_alloy_rho_r = new double*[atom_types_nums];
    EAM_alloy_cubic_spline_coefficient_rho_r = new double*[atom_types_nums];
    total_pair_nums = 0.5 * atom_types_nums * (atom_types_nums + 1);
    EAM_alloy_rfhi_r = new double*[total_pair_nums];
    EAM_alloy_cubic_spline_coefficient_rfhi_r = new double*[total_pair_nums];

    for (int i = 0; i < atom_types_nums; ++i)
    {
        EAM_alloy_F_rho[i] = new double[EAM_Nrho];
        EAM_alloy_rho_r[i] = new double[EAM_Nr];
    }//initialise the array

    for (int i = 0; i < total_pair_nums; ++i)
    {
        EAM_alloy_rfhi_r[i] = new double[EAM_Nr];
    }


    for (int i = 0; i < atom_types_nums; ++i)
    {
        double mass;
        file >> mass;
        file >> mass;// skip the first index

        const double EPSILON = 1e-5;  // you can adjust the number due to your requierment

        if (random)
        {
            if(i == type)
            {
                int nums = atom_0.total_num * percentage;
                std::vector<int> random_index(atom_0.total_num);
                for (int j = 0; j < atom_0.total_num; ++j)
                random_index[j] = j;

                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::default_random_engine rng(seed);
                std::shuffle(random_index.begin(), random_index.end(), rng);// randomise the vector


                for (int j = 0; j < nums; ++j)
                {
                    atom_0.mass[random_index[j]] = mass;
                    atom_0.mass_inv[random_index[j]] = 1 / mass;
                }// randomise the atoms
            }
        }

        for (int n = 0; n < atom_0.total_num; ++n)
        {
            if (fabs(atom_0.mass[n] - mass) < EPSILON)
            atom_0.type[n] = i;
        } // label the atoms



        std::getline(file, line);
        for (int j = 0; j < EAM_Nrho; ++j)
        {
            file >> EAM_alloy_F_rho[i][j];
        }//read the data about Frho

        for (int j = 0; j < EAM_Nr; ++j)
        {
            file >> EAM_alloy_rho_r[i][j];
        }// read the data about rhor

    }   

    for (int i = 0; i < total_pair_nums; ++i)
    {
        for (int j = 0; j < EAM_Nr; ++j)
        {
            file >> EAM_alloy_rfhi_r[i][j];
        }
    }


    file.close();

    EAM_dr_inv = 1 / EAM_dr;
    EAM_drho_inv = 1 / EAM_drho;
    coefficient::cutoff_box_x = coefficient::cutoff_box_y = coefficient::cutoff_box_z = coefficient::cutoff + 1;
    coefficient::cut_off_nl = coefficient::cutoff + coefficient::nl_delta;
    coefficient::cut_off_nl_square = coefficient::cut_off_nl * coefficient::cut_off_nl;
    coefficient::cutoffSquare = coefficient::cutoff * coefficient::cutoff;
}






void EAM_cubic_spline_interpolation_partial(int& total_num, double*& values, double& step, double*& obj)
{
            int num = 4 * (total_num - 1);
            Eigen::SparseMatrix<double> matrix_for_r(num, num);//   the parameter is a, b, c, d with the form as a+bx+cx2+dx3
            Eigen::VectorXd vector_for_r(num);      // the derivative is b + 2cx + 3dx2, 2c + 6dx
            vector_for_r.setZero();

            matrix_for_r.insert(0, 2) = 2;
            matrix_for_r.insert(1, num - 1) = 6 * (total_num - 1) * step;
            matrix_for_r.insert(1, num - 2) = 2;// boundry condition
            matrix_for_r.insert(2, 0) = 1;
            vector_for_r(2) = values[0];
            double r_max = step * (total_num - 1);
            matrix_for_r.insert(3, num - 1) = r_max * r_max * r_max;
            matrix_for_r.insert(3, num - 2) = r_max * r_max;
            matrix_for_r.insert(3, num - 3) = r_max;
            matrix_for_r.insert(3, num - 4) = 1;
            vector_for_r(3) = values[total_num - 1];

            #pragma omp parallel for
            for (int i = 1; i < total_num - 1; ++i)
            {
                double r = step * i;
                int row = i * 4;
                int col = (i - 1) * 4;
                matrix_for_r.insert(row, col) = 1;
                matrix_for_r.insert(row, col + 1) = r;
                matrix_for_r.insert(row, col + 2) = r * r;
                matrix_for_r.insert(row, col + 3) = r * r * r;
                vector_for_r(row) = values[i];  // the first equation

                row += 1;
                matrix_for_r.insert(row, col + 4) = 1;
                matrix_for_r.insert(row, col + 5) = r;
                matrix_for_r.insert(row, col + 6) = r * r;
                matrix_for_r.insert(row, col + 7) = r * r * r;
                vector_for_r(row) = values[i];  // the second equation

                row += 1;
                matrix_for_r.insert(row, col + 1) = 1;
                matrix_for_r.insert(row, col + 2) = 2 * r;
                matrix_for_r.insert(row, col + 3) = 3 * r * r;


                matrix_for_r.insert(row, col + 5) = - 1;
                matrix_for_r.insert(row, col + 6) = - 2 * r;
                matrix_for_r.insert(row, col + 7) = - 3 * r * r;
                // the third equation about the first derivative


                row += 1;
                matrix_for_r.insert(row, col + 2) = 2;
                matrix_for_r.insert(row, col + 3) = 6 * r;

                matrix_for_r.insert(row, col + 6) = - 2;
                matrix_for_r.insert(row, col + 7) = - 6 * r;
                // the forth equation about the second derivative
                
            }


            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(matrix_for_r);

            Eigen::VectorXd ans = solver.solve(vector_for_r);
            //solve the equations



            obj = new double [num];
            for (int i = 0; i < num; ++i)
            {
                obj[i] = ans(i);
            }

}


void potential_energy::EAM_cubic_spline_interpolation()
{
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            EAM_cubic_spline_interpolation_partial(EAM_Nr, EAM_Z_r, EAM_dr, EAM_cubic_spline_coefficient_r_Zr);
        }
        #pragma omp section
        {
            EAM_cubic_spline_interpolation_partial(EAM_Nr, EAM_rho_r, EAM_dr, EAM_cubic_spline_coefficient_r_rhor);
        }
        #pragma omp section
        {
            EAM_cubic_spline_interpolation_partial(EAM_Nrho, EAM_F_rho, EAM_drho, EAM_cubic_spline_coefficient_rho_Frho);
        }
    }
}


void potential_energy::EAM_alloy_cubic_spline_interpolation()
{
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            for (int i = 0; i < atom_types_nums; ++i)
            {
                EAM_cubic_spline_interpolation_partial(EAM_Nrho, EAM_alloy_F_rho[i], EAM_drho, EAM_alloy_cubic_spline_coefficient_F_rho[i]);
                delete[] EAM_alloy_F_rho[i];
            }
        }
        #pragma omp section
        {
            for (int i = 0; i < atom_types_nums; ++i)
            {
                EAM_cubic_spline_interpolation_partial(EAM_Nr, EAM_alloy_rho_r[i], EAM_dr, EAM_alloy_cubic_spline_coefficient_rho_r[i]);
                delete[] EAM_alloy_rho_r[i];
            }
        }
        #pragma omp section
        {
            for (int i = 0; i < total_pair_nums; ++i)
            {
                EAM_cubic_spline_interpolation_partial(EAM_Nr, EAM_alloy_rfhi_r[i], EAM_dr, EAM_alloy_cubic_spline_coefficient_rfhi_r[i]);
                delete[] EAM_alloy_rfhi_r[i];
            }
        }
    }
}





//just as the name says
double potential_energy::get_the_value_of_spline(double*& coefficient_num, double& x, double& step_inv)//the double* point to an array
{
    int n = x * step_inv;
    int index = n * 4;
    return coefficient_num[index] + coefficient_num[index + 1] * x + coefficient_num[index + 2] * x * x + coefficient_num[index + 3] * x * x * x;
}

//just as the name says
double potential_energy::get_the_derivative_of_the_spline(double*& coefficient_num, double& x, double& step_inv)
{
    int n = x * step_inv;
    int index = n * 4;
    return coefficient_num[index + 1] + coefficient_num[index + 2] * 2 * x + coefficient_num[index + 3] * 3 * x * x;
}


void transform_into_unit_vector(double* arr)
{
    double mold = sqrt(arr[0] * arr[0] + arr[1] * arr[1] + arr[2] * arr[2]);
    double mold_inv = 1 / mold;
    arr[0] *= mold_inv;
    arr[1] *= mold_inv;
    arr[2] *= mold_inv;
}


void potential_energy::EAM_potential_cubicspline(atom& atom_0)
{
    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();

    delete[] atom_0.force;
    atom_0.force = new three_dim_vector[atom_0.total_num]();
    //#pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int total_num_plus_1 = atom_0.total_num + 1;
    int* real_neighbour_list = new int[atom_0.total_num * (total_num_plus_1)]();
    std::vector<std::vector<std::vector<double>>> r_all(
        atom_0.total_num, 
        std::vector<std::vector<double>>(
            atom_0.total_num, 
            std::vector<double>(4, 0.0)
            ));   // the 4 element is mold, x, y, z this used for store all the data

    //#pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        for (int j = 0; j < atom_0.total_num; ++j)
        {
            if (j == i)
            continue;
            double xij = atom_0.coordinate[j].x - atom_0.coordinate[i].x;
            double yij = atom_0.coordinate[j].y - atom_0.coordinate[i].y;
            double zij = atom_0.coordinate[j].z - atom_0.coordinate[i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            double old_r[3] = {xij, yij, zij};
            double new_r[3];
            compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
            xij = new_r[0];
            yij = new_r[1];
            zij = new_r[2];
            double rij = sqrt(xij * xij + yij * yij + zij * zij);
            if (rij > coefficient::cutoff)
            continue;
            double rij_inv = 1 / rij;
            xij *= rij_inv;
            yij *= rij_inv;
            zij *= rij_inv;
            int prefix = i * total_num_plus_1;
            real_neighbour_list[prefix] += 1;
            real_neighbour_list[prefix + real_neighbour_list[prefix]] = j;
            r_all[i][j][0] = rij;
            r_all[i][j][1] = xij;
            r_all[i][j][2] = yij;
            r_all[i][j][3] = zij;
            atom_0.rho[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_r_rhor, rij, EAM_dr_inv);
            double Z_r = get_the_value_of_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            atom_0.pe[i] += 27.2 * 0.529 * Z_r * Z_r * rij_inv * 0.5;
            double dZ_r = get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            double force = 27.2 * 0.529 * (2 * Z_r * rij  * dZ_r - Z_r * Z_r) * rij_inv * rij_inv;  // this is the two body force part

            double Fx, Fy, Fz;
            Fx = force * xij;
            Fy = force * yij;
            Fz = force * zij;
            atom_0.force[i].x += Fx;
            atom_0.force[i].y += Fy;
            atom_0.force[i].z += Fz;

            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = rij * xij;
                ry = rij * yij;
                rz = rij * zij;
                Fx *= 0.5;
                Fy *= 0.5;
                Fz *= 0.5;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }

        }
    }

    //#pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * (atom_0.total_num + 1);
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {

            
            double force = ((get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv))
            + (get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[real_neighbour_list[prefix + j]], EAM_drho_inv)))
            * get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_rhor, r_all[i][real_neighbour_list[prefix + j]][0], EAM_dr_inv); //this is the multitude body force part

            double Fx, Fy, Fz;
            Fx = force * r_all[i][real_neighbour_list[prefix + j]][1];
            Fy = force * r_all[i][real_neighbour_list[prefix + j]][2];
            Fz = force * r_all[i][real_neighbour_list[prefix + j]][3];

            atom_0.force[i].x += Fx;
            atom_0.force[i].y += Fy;
            atom_0.force[i].z += Fz;    //calculate the many body force

            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][1];
                ry = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][2];
                rz = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][3];
                Fx *= 0.5;
                Fy *= 0.5;
                Fz *= 0.5;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }

        }
    }
    delete[] real_neighbour_list;
}




void potential_energy::EAM_alloy_potential_cubicspline(atom& atom_0)
{
    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();
    
    //#pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.force[n].x = atom_0.force[n].y = atom_0.force[n].z = atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int total_num_plus_1 = atom_0.total_num + 1;
    int* real_neighbour_list = new int[atom_0.total_num * (total_num_plus_1)]();
    std::vector<std::vector<std::vector<double>>> r_all(
        atom_0.total_num, 
        std::vector<std::vector<double>>(
            atom_0.total_num, 
            std::vector<double>(4, 0.0)
            ));   // the 4 element is mold, x, y, z this used for store all the data

    //#pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        for (int j = 0; j < atom_0.total_num; ++j)
        {
            if (j == i)
            continue;
            double xij = atom_0.coordinate[j].x - atom_0.coordinate[i].x;
            double yij = atom_0.coordinate[j].y - atom_0.coordinate[i].y;
            double zij = atom_0.coordinate[j].z - atom_0.coordinate[i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            double old_r[3] = {xij, yij, zij};
            double new_r[3];
            compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
            xij = new_r[0];
            yij = new_r[1];
            zij = new_r[2];
            double rij = sqrt(xij * xij + yij * yij + zij * zij);
            if (rij > coefficient::cutoff)
            continue;
            double rij_inv = 1 / rij;
            xij *= rij_inv;
            yij *= rij_inv;
            zij *= rij_inv;
            int prefix = i * total_num_plus_1;
            real_neighbour_list[prefix] += 1;
            real_neighbour_list[prefix + real_neighbour_list[prefix]] = j;
            r_all[i][j][0] = rij;
            r_all[i][j][1] = xij;
            r_all[i][j][2] = yij;
            r_all[i][j][3] = zij;
            int bigger, smaller;
            if (atom_0.type[i] > atom_0.type[j])
            {
                bigger = atom_0.type[i];
                smaller = atom_0.type[j];
            }
            else
            {
                bigger = atom_0.type[j];
                smaller = atom_0.type[i];
            }
            int index = 0.5 * bigger * (bigger + 1) + smaller;
            atom_0.rho[i] += get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_rho_r[atom_0.type[j]], rij, EAM_dr_inv);
            double fi_r = get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_rfhi_r[index], rij, EAM_dr_inv);
            atom_0.pe[i] += 0.5 * fi_r * rij_inv;
            double dfi_r = get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rfhi_r[index], rij, EAM_dr_inv);
            double force = (dfi_r - fi_r * rij_inv) * rij_inv;  // this is the two body force part

            double Fx, Fy, Fz;
            Fx = force * xij;
            Fy = force * yij;
            Fz = force * zij;

            atom_0.force[i].x += Fx;
            atom_0.force[i].y += Fy;
            atom_0.force[i].z += Fz;

            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = rij * xij;
                ry = rij * yij;
                rz = rij * zij;
                Fx *= 0.5;
                Fy *= 0.5;
                Fz *= 0.5;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }
        }
    }

    //#pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_F_rho[atom_0.type[i]], atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * (atom_0.total_num + 1);
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {         
            double force = (get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[atom_0.type[i]], atom_0.rho[i], EAM_drho_inv)
            * get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rho_r[atom_0.type[real_neighbour_list[prefix + j]]], r_all[i][real_neighbour_list[prefix + j]][0], EAM_dr_inv))
            + (get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[atom_0.type[real_neighbour_list[prefix + j]]], atom_0.rho[real_neighbour_list[prefix + j]], EAM_drho_inv)
            * get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rho_r[atom_0.type[i]], r_all[i][real_neighbour_list[prefix + j]][0], EAM_dr_inv)); //this is the multitude body force part


            double Fx, Fy, Fz;
            Fx = force * r_all[i][real_neighbour_list[prefix + j]][1];
            Fy = force * r_all[i][real_neighbour_list[prefix + j]][2];
            Fz = force * r_all[i][real_neighbour_list[prefix + j]][3];


            atom_0.force[i].x += Fx;
            atom_0.force[i].y += Fy;
            atom_0.force[i].z += Fz;    //calculate the many body force

            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][1];
                ry = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][2];
                rz = r_all[i][real_neighbour_list[prefix + j]][0] * r_all[i][real_neighbour_list[prefix + j]][3];
                Fx *= 0.5;
                Fy *= 0.5;
                Fz *= 0.5;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }
        }
    }
    delete[] real_neighbour_list;
}






void potential_energy::EAM_nl_linear_interpolation(atom& atom_0)
{
    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();
    //#pragma omp parallel for

    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.force[n].x = atom_0.force[n].y = atom_0.force[n].z = atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int* real_neighbour_list = new int[atom_0.total_num * (coefficient::max_nums_in_cutoff)]();
    vec* r_all = new vec[atom_0.total_num * coefficient::max_nums_in_cutoff];  // the 4 element is mold, x, y, z, this used for store all the data



    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        int prefix = i * coefficient::max_num_in_nl;
        int prefix_of_real_nl = i * coefficient::max_nums_in_cutoff;
        int nums_of_neighbour = atom_0.neighbour_list[prefix];
        for (int j = 1; j <= nums_of_neighbour; ++j)
        {
            int real_j = atom_0.neighbour_list[j + prefix];//get the real index
            double xij = atom_0.coordinate[real_j].x - atom_0.coordinate[i].x;
            double yij = atom_0.coordinate[real_j].y - atom_0.coordinate[i].y;
            double zij = atom_0.coordinate[real_j].z - atom_0.coordinate[i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            if (!atom_0.if_othogonal)
            {
                double old_r[3] = {xij, yij, zij};
                double new_r[3];
                compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
                xij = new_r[0];
                yij = new_r[1];
                zij = new_r[2];
            }
            double rij = sqrt(xij * xij + yij * yij + zij * zij);
            if (rij > coefficient::cutoff)
            continue;
            double rij_inv = 1 / rij;
            xij *= rij_inv;
            yij *= rij_inv;
            zij *= rij_inv;
            real_neighbour_list[prefix_of_real_nl] += 1;
            int put_place = real_neighbour_list[prefix_of_real_nl] + prefix_of_real_nl;//  to align the numbers
            real_neighbour_list[put_place] = real_j;
            r_all[put_place].mold = rij;
            r_all[put_place].x = xij;
            r_all[put_place].y = yij;
            r_all[put_place].z = zij;
            double rho = get_the_value_of_spline(EAM_cubic_spline_coefficient_r_rhor, rij, EAM_dr_inv);
            #pragma omp atomic
            atom_0.rho[real_j] += rho;
            #pragma omp atomic
            atom_0.rho[i] += rho;
            double Z_r = get_the_value_of_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            double pe = 27.2 * 0.529 * Z_r * Z_r * rij_inv * 0.5;
            #pragma omp atomic
            atom_0.pe[i] += pe;
            #pragma omp atomic
            atom_0.pe[real_j] += pe;
            double dZ_r = get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            double force = 27.2 * 0.529 * (2 * Z_r * rij  * dZ_r - Z_r * Z_r) * rij_inv * rij_inv;  // this is the two body force part
            r_all[put_place].force = force;
        }
    }




    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * coefficient::max_nums_in_cutoff;
        double local_force_x = 0.0, local_force_y = 0.0, local_force_z = 0.0;
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {
            int index_j = real_neighbour_list[prefix + j];
            double force = ((get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv))
            + (get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[index_j], EAM_drho_inv)))
            * get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_rhor, r_all[prefix + j].mold, EAM_dr_inv)
            + r_all[prefix + j].force; //this is the multitude body force part

            double Fx, Fy, Fz;
            Fx = force * r_all[prefix + j].x;
            Fy = force * r_all[prefix + j].y;
            Fz = force * r_all[prefix + j].z;

            local_force_x += Fx;
            local_force_y += Fy;
            local_force_z += Fz;

            #pragma omp atomic
            atom_0.force[index_j].x -= Fx;
        
            #pragma omp atomic
            atom_0.force[index_j].y -= Fy;
        
            #pragma omp atomic
            atom_0.force[index_j].z -= Fz; //calculate the many body force


            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = r_all[prefix + j].mold * r_all[prefix + j].x;
                ry = r_all[prefix + j].mold * r_all[prefix + j].y;
                rz = r_all[prefix + j].mold * r_all[prefix + j].z;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }

        }
        #pragma omp atomic
        atom_0.force[i].x += local_force_x;

        #pragma omp atomic
        atom_0.force[i].y += local_force_y;

        #pragma omp atomic
        atom_0.force[i].z += local_force_z;

    }
    delete[] real_neighbour_list;
    delete[] r_all;
}


void potential_energy::EAM_alloy_nl_cubic_spline(atom& atom_0)
{
    if (atom_0.virial_switch)
    atom_0.clean_the_virial_tensor();

    //#pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.force[n].x = atom_0.force[n].y = atom_0.force[n].z = atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int* real_neighbour_list = new int[atom_0.total_num * (coefficient::max_nums_in_cutoff)]();
    vec* r_all = new vec[atom_0.total_num * coefficient::max_nums_in_cutoff];  // the 4 element is mold, x, y, z, this used for store all the data



    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        int prefix = i * coefficient::max_num_in_nl;
        int prefix_of_real_nl = i * coefficient::max_nums_in_cutoff;
        int nums_of_neighbour = atom_0.neighbour_list[prefix];
        int type_i = atom_0.type[i];
        for (int j = 1; j <= nums_of_neighbour; ++j)
        {
            int real_j = atom_0.neighbour_list[j + prefix];//get the real index
            double xij = atom_0.coordinate[real_j].x - atom_0.coordinate[i].x;
            double yij = atom_0.coordinate[real_j].y - atom_0.coordinate[i].y;
            double zij = atom_0.coordinate[real_j].z - atom_0.coordinate[i].z;
            mirror_constrain(atom_0.box, xij, yij, zij);
            if (!atom_0.if_othogonal)
            {
                double old_r[3] = {xij, yij, zij};
                double new_r[3];
                compute_1dv_times_3dm(old_r, atom_0.basic_vectors, new_r);
                xij = new_r[0];
                yij = new_r[1];
                zij = new_r[2];
            }
            double rij = sqrt(xij * xij + yij * yij + zij * zij);

            if (rij > coefficient::cutoff)
            continue;
            double rij_inv = 1 / rij;
            xij *= rij_inv;
            yij *= rij_inv;
            zij *= rij_inv;
            real_neighbour_list[prefix_of_real_nl] += 1;
            int put_place = real_neighbour_list[prefix_of_real_nl] + prefix_of_real_nl;//  to align the numbers
            real_neighbour_list[put_place] = real_j;
            r_all[put_place].mold = rij;
            r_all[put_place].x = xij;
            r_all[put_place].y = yij;
            r_all[put_place].z = zij;


            int type_j = atom_0.type[real_j];   //get the type


            int index;
            if (type_i > type_j)
            {
                index = 0.5 * type_i * (type_i + 1) + type_j;
            }
            else
            {
                index = 0.5 * type_j * (type_j + 1) + type_i;
            }


            double rho = get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_rho_r[type_j], rij, EAM_dr_inv);
            #pragma omp atomic
            atom_0.rho[i] += rho;
            if (type_i != type_j)
            {
                rho = get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_rho_r[type_i], rij, EAM_dr_inv);
            }
            #pragma omp atomic
            atom_0.rho[real_j] += rho;  //add the rho

            double fi_r = get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_rfhi_r[index], rij, EAM_dr_inv);
            double pe = 0.5 * fi_r * rij_inv;
            #pragma omp atomic
            atom_0.pe[i] += pe;
            #pragma omp atomic
            atom_0.pe[real_j] += pe;
            double dfi_r = get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rfhi_r[index], rij, EAM_dr_inv);
            r_all[put_place].force = (dfi_r - fi_r * rij_inv) * rij_inv;   // this is the two body force part

        }
    }


    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_alloy_cubic_spline_coefficient_F_rho[atom_0.type[i]], atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * coefficient::max_nums_in_cutoff;
        int type_i = atom_0.type[i];
        double local_force_x = 0.0, local_force_y = 0.0, local_force_z = 0.0;
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {
            int index_j = real_neighbour_list[prefix + j];
            double force;
            int type_j = atom_0.type[index_j];   //get the type
            if (type_i != type_j)
            {
                force = (get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[type_i], atom_0.rho[i], EAM_drho_inv)
                * get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rho_r[type_j], r_all[prefix + j].mold, EAM_dr_inv))
                + (get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[type_j], atom_0.rho[index_j], EAM_drho_inv)
                * get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rho_r[type_i], r_all[prefix + j].mold, EAM_dr_inv))
                + r_all[prefix + j].force; //this is the multitude body force part
            }
            else
            {
                force = ((get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[type_i], atom_0.rho[i], EAM_drho_inv)
                + get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_F_rho[type_j], atom_0.rho[index_j], EAM_drho_inv))
                * get_the_derivative_of_the_spline(EAM_alloy_cubic_spline_coefficient_rho_r[type_i], r_all[prefix + j].mold, EAM_dr_inv))
                + r_all[prefix + j].force; //this is the multitude body force part
            }

            double Fx, Fy, Fz;
            Fx = force * r_all[prefix + j].x;
            Fy = force * r_all[prefix + j].y;
            Fz = force * r_all[prefix + j].z;



            local_force_x += Fx;
            local_force_y += Fy;
            local_force_z += Fz;

            #pragma omp atomic
            atom_0.force[index_j].x -= Fx;
        
            #pragma omp atomic
            atom_0.force[index_j].y -= Fy;
        
            #pragma omp atomic
            atom_0.force[index_j].z -= Fz;//calculate the many body force

            if (atom_0.virial_switch)
            {
                double rx, ry, rz;
                rx = r_all[prefix + j].mold * r_all[prefix + j].x;
                ry = r_all[prefix + j].mold * r_all[prefix + j].y;
                rz = r_all[prefix + j].mold * r_all[prefix + j].z;
                atom_0.update_the_virial_tensor(Fx, Fy, Fz, rx, ry, rz, i);
            }
        }

        #pragma omp atomic
        atom_0.force[i].x += local_force_x;

        #pragma omp atomic
        atom_0.force[i].y += local_force_y;

        #pragma omp atomic
        atom_0.force[i].z += local_force_z;
    }

    delete[] real_neighbour_list;
    delete[] r_all;
}



