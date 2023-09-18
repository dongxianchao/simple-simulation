#include "../../include/atom/atom_potential.h"
#include <iostream>


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



atom::atom(std::vector<int> numcells, int numbers,
    std::vector<double> mass_0,
    std::vector<double> x_0,
    std::vector<double> y_0,
    std::vector<double> z_0, cell& cell_0)
    {
        cell_num = numbers;
        total_num = cell_num * numcells[0] * numcells[1] * numcells[2];
        mass.resize(total_num, 0);
        mass_inv.resize(total_num, 0);
        x.resize(total_num, 0);
        y.resize(total_num, 0);
        z.resize(total_num, 0);
        vx.resize(total_num, 0);
        vy.resize(total_num, 0);
        vz.resize(total_num, 0);
        fx.resize(total_num, 0);
        fy.resize(total_num, 0);
        fz.resize(total_num, 0);
        pe.resize(total_num, 0);
        x0.resize(total_num, 0);
        y0.resize(total_num, 0);
        z0.resize(total_num, 0); // allocate memory
        rho.resize(total_num, 0);




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
                        x[n] = (ix + x_0[i]) * cell_0.lattice_x;
                        y[n] = (iy + y_0[i]) * cell_0.lattice_y;
                        z[n] = (iz + z_0[i]) * cell_0.lattice_z;
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
            energy += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) * mass[i];
        }
        kineticenergy = energy;
        return energy;
    }

    double atom::find_temperature(double& energy)
    {
        double t = 2 * energy / (3 * coefficient::k_B * total_num);
        return t;
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
            vx[n] = -1 + (rand() * 2.0) / RAND_MAX;
            vy[n] = -1 + (rand() * 2.0) / RAND_MAX;
            vz[n] = -1 + (rand() * 2.0) / RAND_MAX;
            centerofmassvelocity[0] += vx[n] * mass[n];
            centerofmassvelocity[1] += vy[n] * mass[n];
            centerofmassvelocity[2] += vz[n] * mass[n];
        }
        centerofmassvelocity[0] /= total_mass;
        centerofmassvelocity[1] /= total_mass;
        centerofmassvelocity[2] /= total_mass;



        for (int i = 0; i < total_num; ++i) //compute the relative velocity
        {
            vx[i] -= centerofmassvelocity[0];
            vy[i] -= centerofmassvelocity[1];
            vz[i] -= centerofmassvelocity[2];
        }



        double energy = kinetic_energy();// calculate the temperature
        double t = find_temperature(energy);

        double scalar = sqrt(T0 / t);       //find the velocity scalar



        for (int i = 0; i < total_num; ++i)     //update the velocity
        {
            vx[i] *= scalar;
            vy[i] *= scalar;
            vz[i] *= scalar;
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
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.fx[n] = atom_0.fy[n] = atom_0.fz[n] = atom_0.pe[n] = 0;
    }

    for (int i = 0; i < atom_0.total_num - 1; ++i)  //calculate the force
    {
        for (int j = i + 1; j < atom_0.total_num; ++j)
        {
            double xij = atom_0.x[j] - atom_0.x[i];
            double yij = atom_0.y[j] - atom_0.y[i];
            double zij = atom_0.z[j] - atom_0.z[i];
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
            atom_0.fx[i] += f_ij * xij;
            atom_0.fx[j] -= f_ij * xij;
            atom_0.fy[i] += f_ij * yij;
            atom_0.fy[j] -= f_ij * yij;
            atom_0.fz[i] += f_ij * zij;
            atom_0.fz[j] -= f_ij * zij;
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
    return kineticenergy + potentialenergy;
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


    coefficient::max_num_in_boxes = (coefficient::cutoff_box_x / cell_0.lattice_x) *
    (coefficient::cutoff_box_y / cell_0.lattice_y) * (coefficient::cutoff_box_z / cell_0.lattice_z)
    * (cell_num * 16);
    coefficient::max_num_in_nl = 4 * coefficient::max_num_in_boxes;
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
    boxes[prefix] += 1;
    boxes[boxes[prefix] + prefix] = index;
}


void atom::insert_in_boxes(int (&place)[3], int index)
{
    int prefix = (place[0] * boxes_y * boxes_z + place[1] * boxes_z + place[2])
    * coefficient::max_num_in_boxes;
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
    #pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.fx[n] = atom_0.fy[n] = atom_0.fz[n] = atom_0.pe[n] = 0;
    }// initialise the force and potential energy

    #pragma omp parallel for
    for (int index_i = 0; index_i < atom_0.total_num; ++index_i)  //calculate the force
    {
        int prefix = index_i * coefficient::max_num_in_nl;
        for(int s = 1; s <= atom_0.neighbour_list[prefix]; ++s)
        {
            int index_j = atom_0.neighbour_list[prefix + s];
            double xij = atom_0.x[index_j] - atom_0.x[index_i];
            double yij = atom_0.y[index_j] - atom_0.y[index_i];
            double zij = atom_0.z[index_j] - atom_0.z[index_i];
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
            atom_0.pe[index_i] += 0.5 * (coefficient::e4s12 * r12inv - coefficient::e4s6 * r6inv);
            atom_0.fx[index_i] += f_ij * xij;
            atom_0.fy[index_i] += f_ij * yij;
            atom_0.fz[index_i] += f_ij * zij;
        }
    }
}



void atom::initialise_neighbour_list()// used to initialise nl
{
    neighbour_list = new int[coefficient::max_num_in_nl * total_num]();
}

void atom::update_neighbour_list()
{
    #pragma omp parallel for
    for (int index_i = 0; index_i < total_num; ++index_i)  //   for every atom
    {
        int place[3];
        place[0] = x[index_i] * coefficient::cutoff_box_inv_x;
        place[1] = y[index_i] * coefficient::cutoff_box_inv_y;
        place[2] = z[index_i] * coefficient::cutoff_box_inv_z;

        if (place[0] == boxes_x) place[0] -= 1;
        if (place[1] == boxes_y) place[1] -= 1;
        if (place[2] == boxes_z) place[2] -= 1;//    the same function as the get_place function

        for (int i = place[0] - 1; i <= place[0] + 1; ++i)
        {
            for (int j = place[1] - 1; j <= place[1] + 1; ++j)
            {
                for (int k = place[2] - 1; k <= place[2] + 1; ++k)
                {
                    int x_index = i;
                    int y_index = j;
                    int z_index = k;
                    if (x_index == -1 || x_index == boxes_x)
                    x_index = (x_index == -1)? boxes_x -1 : 0;
                    if (y_index == -1 || y_index == boxes_y)
                    y_index = (y_index == -1)? boxes_y -1 : 0;
                    if (z_index == -1 || z_index == boxes_z)
                    z_index = (z_index == -1)? boxes_z -1 : 0;

                    int prefix = (x_index * boxes_y * boxes_z
                     + y_index * boxes_z + z_index) * coefficient::max_num_in_boxes;
                    int num = boxes[prefix];
                    for (int s = 1; s <= num; ++s)
                    {
                        int index_j = boxes[prefix + s];
                        if (index_j == index_i)
                        continue;// prevent the same atom
                        double xij = x[index_j] - x[index_i];
                        double yij = y[index_j] - y[index_i];
                        double zij = z[index_j] - z[index_i];
                        mirror_constrain(box, xij, yij, zij);
                        double old_r[3] = {xij, yij, zij};
                        double new_r[3];
                        compute_1dv_times_3dm(old_r, basic_vectors, new_r);
                        xij = new_r[0];
                        yij = new_r[1];
                        zij = new_r[2];
                        double r2 = xij*xij + yij*yij + zij*zij;

                        if (r2 > coefficient::cut_off_nl_square)
                        continue;// set the cutoff

                        int pre = coefficient::max_num_in_nl * index_i;
                        neighbour_list[pre] += 1;
                        neighbour_list[neighbour_list[pre] + pre] = index_j;
                    }
                }
            }
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

void potential_energy::get_the_EAM_data(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    file >> EAM_Nrho >> EAM_drho >> EAM_Nr >> EAM_dr >> coefficient::cutoff;


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
            vector_for_r(3) = values[num - 1];

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
    #pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.fx[n] = atom_0.fy[n] = atom_0.fz[n] = atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int total_num_plus_1 = atom_0.total_num + 1;
    int* real_neighbour_list = new int[atom_0.total_num * (total_num_plus_1)]();
    std::vector<std::vector<std::vector<double>>> r_all(
        atom_0.total_num, 
        std::vector<std::vector<double>>(
            atom_0.total_num, 
            std::vector<double>(4, 0.0)
            ));   // the 4 element is mold, x, y, z, this used for store all the data

    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        for (int j = 0; j < atom_0.total_num; ++j)
        {
            if (j == i)
            continue;
            double xij = atom_0.x[j] - atom_0.x[i];
            double yij = atom_0.y[j] - atom_0.y[i];
            double zij = atom_0.z[j] - atom_0.z[i];
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
            double force = 2 * 27.2 * 0.529 * (2 * Z_r * rij  * dZ_r - Z_r * Z_r) * rij_inv * rij_inv;  // this is the two body force part
            atom_0.fx[i] += force * xij;
            atom_0.fy[i] += force * yij;
            atom_0.fz[i] += force * zij;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * (atom_0.total_num + 1);
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {

            
            double force = ((get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv))
            + (get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[real_neighbour_list[prefix + j]], EAM_drho_inv)))
            * get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_rhor, r_all[i][real_neighbour_list[prefix + j]][0], EAM_dr_inv); //this is the multitude body force part


            atom_0.fx[i] += force * r_all[i][real_neighbour_list[prefix + j]][1];
            atom_0.fy[i] += force * r_all[i][real_neighbour_list[prefix + j]][2];
            atom_0.fz[i] += force * r_all[i][real_neighbour_list[prefix + j]][3];    //calculate the many body force

        }
    }
    delete[] real_neighbour_list;
}




void potential_energy::EAM_nl_linear_interpolation(atom& atom_0)
{
    #pragma omp parallel for
    for (int n = 0; n < atom_0.total_num; ++n)
    {
        atom_0.fx[n] = atom_0.fy[n] = atom_0.fz[n] = atom_0.pe[n] = atom_0.rho[n] = 0;
    }// initialise the data


    int* real_neighbour_list = new int[atom_0.total_num * (coefficient::max_num_in_nl)]();
    std::vector<std::vector<std::vector<double>>> r_all(
        atom_0.total_num, 
        std::vector<std::vector<double>>(
            coefficient::max_num_in_nl, 
            std::vector<double>(4, 0.0)
            ));   // the 4 element is mold, x, y, z, this used for store all the data

    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        int prefix = i * coefficient::max_num_in_nl;
        int nums_of_neighbour = atom_0.neighbour_list[prefix];
        for (int j = 1; j <= nums_of_neighbour; ++j)
        {
            double real_j = atom_0.neighbour_list[j + prefix];//get the real index
            double xij = atom_0.x[real_j] - atom_0.x[i];
            double yij = atom_0.y[real_j] - atom_0.y[i];
            double zij = atom_0.z[real_j] - atom_0.z[i];
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
            real_neighbour_list[prefix] += 1;
            int put_place = real_neighbour_list[prefix];//  to align the numbers
            real_neighbour_list[prefix + put_place] = real_j;
            r_all[i][put_place][0] = rij;
            r_all[i][put_place][1] = xij;
            r_all[i][put_place][2] = yij;
            r_all[i][put_place][3] = zij;
            atom_0.rho[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_r_rhor, rij, EAM_dr_inv);
            double Z_r = get_the_value_of_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            atom_0.pe[i] += 27.2 * 0.529 * Z_r * Z_r * rij_inv * 0.5;
            double dZ_r = get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_Zr, rij, EAM_dr_inv);
            double force = 2 * 27.2 * 0.529 * (2 * Z_r * rij  * dZ_r - Z_r * Z_r) * rij_inv * rij_inv;  // this is the two body force part
            atom_0.fx[i] += force * xij;
            atom_0.fy[i] += force * yij;
            atom_0.fz[i] += force * zij;
        }
    }


    #pragma omp parallel for
    for (int i = 0; i < atom_0.total_num; ++i)
    {
        atom_0.pe[i] += get_the_value_of_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv);// calculate the potential energy
        int prefix = i * coefficient::max_num_in_nl;
        for (int j = 1; j <= real_neighbour_list[prefix]; ++j)
        {
            double index_j = real_neighbour_list[prefix + j];
            double force = ((get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[i], EAM_drho_inv))
            + (get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_rho_Frho, atom_0.rho[index_j], EAM_drho_inv)))
            * get_the_derivative_of_the_spline(EAM_cubic_spline_coefficient_r_rhor, r_all[i][j][0], EAM_dr_inv); //this is the multitude body force part


            atom_0.fx[i] += force * r_all[i][j][1];
            atom_0.fy[i] += force * r_all[i][j][2];
            atom_0.fz[i] += force * r_all[i][j][3];    //calculate the many body force

        }
    }
    delete[] real_neighbour_list;
}