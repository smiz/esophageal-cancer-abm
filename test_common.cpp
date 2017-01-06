#include "common.h"
#include <cassert>
#include <iostream>
using namespace std;

const double dx = 1.0;
const double stem_cell_density = 20.0;
const int nx = 10, ny = 20, nz = 30;

void test_exponential()
{
	cout << "TEST EXPONENTIAL" << endl;
	double mean = 1.0;
	double sum = 0.0;
	int count = 10000000;
	for (int i = 0; i < count; i++)
	{
		sum += Parameters::getInstance()->exponential(mean);
	}
	double test_mean = (sum/(double)count);
	cout << "Got mean = " << test_mean << endl;
	assert(fabs(test_mean-mean) < 1E-4);
	cout << "TEST PASSED" << endl;
}

void test_direction()
{
	cout << "TEST DIRECTION" << endl;
	int sum = 0;
	int count = 100000000;
	int dx, dy, dz;
	for (int i = 0; i < count; i++)
	{
		Parameters::getInstance()->direction(dx,dy,dz);
		assert(abs(dx)+abs(dy)+abs(dz) == 1);
		sum += dx+dy+dz;
	}
	double mean = (double)(sum)/(double)(count);
	cout << "Mean direction = " << mean << endl;
	assert(fabs(mean) < 1E-3);
	cout << "TEST PASSED" << endl;
}

void test_mutation_interval()
{
	cout << "TEST MUTATION INTERVAL" << endl;
	for (int i = 0; i < NUM_CELL_TYPES; i++)
	{
		assert(Parameters::getInstance()->get_mutation_interval(i) == adevs_inf<double>());
		Parameters::getInstance()->set_mutations_per_year(i+1,i);
		double interval = Parameters::getInstance()->get_mutation_interval(i);
		cout << "type " << i << " interval is " << interval << " years" << endl;
		assert(interval == 1.0/(stem_cell_density*dx*dx*(double)(i+1)));
	}
	cout << "TEST PASSED" << endl;
}

void test_diffusion_rate()
{
	cout << "TEST EXPAND INTERVAL" << endl;
	Parameters::getInstance()->set_diffusion_rate(2.0);
	assert(Parameters::getInstance()->get_expand_interval() == (dx*dx/2.0));
	cout << "TEST PASSED" << endl;
}

void test_wrap()
{
	cout << "TEST WRAP" << endl;
	int x, y, z;
	Parameters* p = Parameters::getInstance();
	x = nx/2; y = ny/2; z = nz+1;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == ny/2 && z == nz+1);
	x = nx/2; y = ny+1; z = nz/2;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == ny+1 && z == nz/2);
	x = nx/2; y = ny/2; z = -1;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == ny/2 && z == -1);
	x = nx/2; y = -1; z = nz/2;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == -1 && z == nz/2);
	x = nx/2; y = ny/2; z = nz/2;
	assert(p->wrap(x,y,z));
	assert(x == nx/2 && y == ny/2 && z == nz/2);
	x = -1; y = ny/2; z = nz/2;
	assert(p->wrap(x,y,z));
	assert(x == nx-1 && y == ny/2 && z == nz/2);
	x = nx; y = ny/2; z = nz/2;
	assert(p->wrap(x,y,z));
	assert(x == 0 && y == ny/2 && z == nz/2);
	x = nx/2; y = ny/2; z = nz;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == ny/2 && z == nz);
	x = nx/2; y = ny; z = nz/2;
	assert(!p->wrap(x,y,z));
	assert(x == nx/2 && y == ny && z == nz/2);
	cout << "TEST PASSED" << endl;
}

int main()
{
	Parameters* p = Parameters::getInstance();
	assert(p != NULL);
	p->set_seed(100);
	p->cell_size(dx);
	p->set_stem_cells_per_mm2(stem_cell_density);
	p->xdim(nx);
	p->ydim(ny);
	p->zdim(nz);
	assert(p->xdim() == nx);
	assert(p->ydim() == ny);
	assert(p->zdim() == nz);
	assert(p->cell_size() == dx);
	test_exponential();
	test_direction();
	test_mutation_interval();
	test_diffusion_rate();
	test_wrap();
	return 0;
}
