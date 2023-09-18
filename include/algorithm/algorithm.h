#include "../atom/atom_potential.h"
class algorithm
{
private:

static inline void velvet_for_update_v(atom& atom_0, double& t); 
//auxiliary for velvet to compute velocity

static inline void velvet_for_update_place(double & dt, atom& atom_0);
//auxiliary for velvet to compute place

static inline void velvet_for_update_place_with_nl(double & dt, atom& atom_0);

public:

static void velvet_LJ(double & dt, atom &atom_0);
static void velvet_LJ_nl(double & dt, atom &atom_0);
static void velvet_EAM(double & dt, atom &atom_0);
static void velvet_EAM_nl(double & dt, atom &atom_0);

};
