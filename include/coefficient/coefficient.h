
class coefficient
{
    //inside the machine, the unit of the energy is 1.6e-19, lemgth is 1e-10, mass is 1.66e-27
    //and temperature is K.
public:
static constexpr double e_energy = 1.6e-19;
static constexpr double angstrom = 1e-10;
static constexpr double amu = 1.66e-27;
static constexpr double k_B = 8.617343e-5;
static constexpr double unit_t_fs = 1.018051e+01;
static constexpr double unit_t_ps = 1.018051e-02;
static constexpr double unit_p_bar = 1602176.565;
static double e48s12;
static double e4s12;
static double e24s6;
static double e4s6;
static double cutoff;
static double cutoffSquare;
static double cutoff_box_x;
static double cutoff_box_y;
static double cutoff_box_z;
static double cutoff_box_inv_x;
static double cutoff_box_inv_y;
static double cutoff_box_inv_z;
static int max_nums_in_cutoff;  //used for determine the max numbers in the neighbour list
static int max_num_in_boxes;
static int max_num_in_nl;
static double nl_delta;
static double nl_delta_square;
static double cut_off_nl;
static double cut_off_nl_square;

static void get_LJ_coefficient(double epsilon,double sigma,double cutoff);
};
