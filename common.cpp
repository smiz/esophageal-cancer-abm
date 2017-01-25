#include "common.h"
#include <fstream>
#include <string>
#include <sstream>

Parameters* Parameters::inst = NULL; // The singleton instance

double Parameters::uniform()
{
	return gsl_rng_uniform(r);
}

double Parameters::exponential(double mu)
{
	return gsl_ran_exponential(r,mu);
}

double Parameters::normal(double mean, double std_dev)
{
	return gsl_ran_gaussian(r,sqrt(std_dev))+mean;
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

void Parameters::direction(int& dx, int& dy)
{
	int dir = gsl_rng_uniform_int (r,4);
	dx = dy = 0;
	if (dir == 0) dx = 1;
	else if (dir == 1) dx = -1;
	else if (dir == 2) dy = 1;
	else dy = -1;
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
	nx(-1),
	ny(-1),
	nz(-1),
	dx(-1.0),
	diff_time(adevs_inf<double>()),
	stem_cells_per_mm2(-1.0),
	be_onset(-1.0)
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
	double mutate_normal = -1.0,
		mutate_be = -1.0,
		mutate_dysplasia = -1.0;
	ifstream fin(filename);
	if (fin.bad())
	{
		cout << "Could not open file " << filename << endl;
		exit(0);
	}
	for (string line; getline(fin,line); !fin.eof())
	{
		istringstream sin(line);
		string param; double value;
		sin >> param >> value;
		if (param == "diffusion_rate")
			set_diffusion_rate(value);	
		else if (param == "mutate_normal")
			mutate_normal = value;
		else if (param == "mutate_be")
			mutate_be = value;
		else if (param == "mutate_dysplasia")
			mutate_dysplasia = value;
		else if (param == "stem_cell_density")
			set_stem_cells_per_mm2(value);
		else if (param == "be_onset_age")
			be_onset_age(value);
		else
		{
			cout << "Unknown parameter " << param << endl;
			exit(0);
		}
	}
	if (be_onset < 0.0)
	{
		cout << "The be_onset_age must be positive." << endl;
		exit(0);
	}
	if (stem_cells_per_mm2 < 0.0)
	{
		cout << "The stem_cell_density must be positive." << endl;
		exit(0);
	}
	if (mutate_normal > 0.0)
		set_mutations_per_year(mutate_normal,NORMAL);
	if (mutate_be > 0.0)
		set_mutations_per_year(mutate_be,BE);
	if (mutate_dysplasia > 0.0)
		set_mutations_per_year(mutate_dysplasia,DYSPLASIA);
	fin.close();
}
