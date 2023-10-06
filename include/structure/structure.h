#include<vector>
#include <cmath>
#include "../../lib/Eigen/Dense"
#include "../../lib/Eigen/Sparse"
#include <unordered_map>
#include <chrono>
#include <map>

class cell // This is used for primitive cell
{
    public:
    std::vector<std::vector<double>> matrix;    // basic vectors
    double lattice_x, lattice_y, lattice_z;     //lattice coefficiency

    cell(double lattice_a, double lattice_b, double lattice_c,
    std::vector<std::vector<double>> mtx = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}})
    : lattice_x(lattice_a), lattice_y(lattice_b), lattice_z(lattice_c){
        for (auto& v : mtx)
        {
        double magnitude = 0.0;
        for (double val : v) {
            magnitude += val * val;
        }
    magnitude = std::sqrt(magnitude);

    v[0] /= magnitude;
    v[1] /= magnitude;
    v[2] /= magnitude;
    }
    matrix = mtx;
    }
};