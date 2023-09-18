#include "../../include/coefficient/coefficient.h"



    double coefficient::e48s12 = 0.0;
    double coefficient::e4s12 = 0.0;
    double coefficient::e24s6 = 0.0;
    double coefficient::e4s6 = 0.0;
    double coefficient::cutoff = 0.0;
    double coefficient::cutoffSquare = 0.0;
    double coefficient::cutoff_box_x = 0.0;
    double coefficient::cutoff_box_inv_x = 0.0;
    double coefficient::cutoff_box_y = 0.0;
    double coefficient::cutoff_box_inv_y = 0.0;
    double coefficient::cutoff_box_z = 0.0;
    double coefficient::cutoff_box_inv_z = 0.0;
    int coefficient::max_num_in_boxes = 0.0;
    int coefficient::max_num_in_nl = 0.0;
    double coefficient::nl_delta = 1.0;
    double coefficient::nl_delta_square = 1.0;
    double coefficient::cut_off_nl = 0.0;
    double coefficient::cut_off_nl_square = 0.0;


void coefficient::get_LJ_coefficient(double epsilon,double sigma,double cutoff_0)
{
    cutoff = cutoff_0;
    cut_off_nl = cutoff + 1;
    cut_off_nl_square = cut_off_nl * cut_off_nl;
    cutoff_box_x= cutoff_box_y = cutoff_box_z = cutoff + 1;
    cutoffSquare = cutoff * cutoff;
    const double sigma3 = sigma * sigma * sigma;
    const double sigma6 = sigma3 * sigma3;
    const double sigma12 = sigma6 * sigma6;
    e24s6 = 24.0 * epsilon * sigma6;
    e48s12 = 48.0 * epsilon * sigma12;
    e4s6 = 4.0 * epsilon * sigma6;
    e4s12 = 4.0 * epsilon * sigma12;
}


