#include "common.h"
#include <fstream>
#include <string>
#include <sstream>

Parameters* Parameters::inst = NULL; // The singleton instance

double Parameters::exponential(double mu)
{
	return gsl_ran_exponential(r,mu);
}

void Parameters::direction(int& dx, int& dy, int& dz)
{
	int dir = gsl_rng_uniform_int (r,6);
	dx = dy = dz = 0;
	if (dir == 0) dx = 1;
	else if (dir == 1) dx = -1;
	else if (dir == 2) dy = 1;
	else if (dir == 3) dy = -1;
	else if (dir == 4) dz = 1;
	else dz = -1;
}

bool Parameters::wrap(int& x, int& y, int& z) const
{
	if (x < 0) x = nx-1;
	else if (x >= nx) x = 0;
	return (y >= 0 && y < ny && z >= 0 && z < nz);
}

void Parameters::set_seed(unsigned long seed)
{
	gsl_rng_set(r,seed);
}

Parameters* Parameters::getInstance()
{
	if (inst == NULL)
		inst = new Parameters();
	return inst;
}

void Parameters::deleteInstance()
{
	delete inst;
	inst = NULL;
}

Parameters::Parameters():
	dx(-1.0),
	stem_cells_per_mm2(-1.0),
	diff_time(adevs_inf<double>()),
	nx(-1),
	ny(-1),
	nz(-1)
{
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,0);
	for (int i = 0; i < NUM_CELL_TYPES; i++)
	{
		mutate_time[i] = adevs_inf<double>();
	}
}

Parameters::~Parameters()
{
	gsl_rng_free(r);
}

void Parameters::load_from_file(const char* filename)
{
	ifstream fin(filename);
	if (fin.bad())
	{
		cout << "Could not open file " << filename << endl;
		exit(0);
	}
	for (string line; getline(fin,line); )
	{
		istringstream sin(line);
		string param; double value;
		sin >> param >> value;
		if (param == "diffusion_rate")
			set_diffusion_rate(value);	
		else if (param == "mutate_normal")
			set_mutations_per_year(value,NORMAL);
		else if (param == "mutate_be")
			set_mutations_per_year(value,BE);
		else if (param == "mutate_dysplasia")
			set_mutations_per_year(value,DYSPLASIA);
		else if (param == "stem_cell_density")
			set_stem_cells_per_mm2(value);
		else
		{
			cout << "Unknown parameter " << param << endl;
			exit(0);
		}
	}
	fin.close();
}
