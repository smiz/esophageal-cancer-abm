#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "common.h"
#include "des/des.h"
#include "des/Cell.h"
#include "des/CellGrid.h"
using namespace std;
using namespace adevs;

/********************************************************************************
Basis Information:
The diameter of your esophagus is 24mm, the diameter of a quarter.
The thickness of your esophageal wall is 4mm, the width of 3 pennies.

We will use 24 mm in our simulation which gives the circumference to be 75.4 mm.  
The grid width is 0.42 mm for ni = 179 or 0.21 mm for ni = 359.
Jumbo biopsy is about 5mmx3mm (12 by 7 grid) or (24 by 14 grid).
********************************************************************************/

#define NUM_LAYERS 5
// Fraction thickness of each layer
static const double LayerThicknessFraction[NUM_LAYERS] =
{
	0.2, // Epithelium
	0.1, // Basement membrane
	0.2, // Lamina propia
	0.1, // Mucular mucosaa
	0.4  // Submucosa
};

static const double grid_size = 0.42;
static const double thickness = 4.0;
static const double circumference = 75.4;
static const double length = circumference;
static const int ni = (circumference / grid_size)+1; // Spatial points in X direction. 
static const int nj = (length / grid_size)+1; // Spatial points in Y direction. 
static const int nk = (thickness / grid_size)+1;	 // Spatial points in Z direction. 

/**
 * This data can be visualized using paraview. See 
 * www.paraview.org/Wiki/ParaView/Data_formats
 */
void PrintCSV(int seq_num, double t)
{
	assert(NUM_CELL_TYPES == 4);
	static const char* names[NUM_CELL_TYPES] = {
		"normal",
		"BE",
		"dysplasia",
		"cancer",
	};
	int types[6] = { 0, 0, 0, 0, 0, 0 };
	char filename[100];
	sprintf(filename,"tumor.csv.%d",seq_num);
	ofstream fout(filename);
	fout << "xcoord,ycoord,zcoord,type" << endl;
	for (int i = 0; i < ni; i++)
		for (int j = 0; j < nj; j++)
			for (int k = 0; k < nk; k++)
			{
				int type = tissue->get_iType(i,j,k);
				types[type]++;
				if (type > 1)
				{
					fout << i << "," << j << "," << k << "," << type << endl;
				}
			}
	fout.close();
	cout << "t = " << t << endl;
	for (int i = 0; i < 6; i++)
	{
		cout << names[i] << " : " << types[i] << " " << endl;
	}
}

//===========================================================================//
bool InitModelCondition(void)
{
	des_model = new des::CellGrid(ni,nj,nk);
	// Set lifespans
	assert(tab_cell_death_rate.NumRows() == 6);
	for (int i = 0; i < 6; i++)
	{
		double rate = tab_cell_death_rate.GetValue(i);
		if (rate > 0.0) rate = DAYS_TO_YEARS/rate;
		des::CellComponent::setCellLifespan(i,rate);
	}
	// Set mitotic period
	assert(tab_mitotic_period.NumRows() == 6);
	for (int i = 0; i < 6; i++)
	{
		des::CellComponent::setMitosisPeriod(i,tab_mitotic_period.GetValue(i)*DAYS_TO_YEARS);
	}
	// Setup mutation data
	assert(tab_transition.NumRows() == 6);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			des::CellComponent::setMutationRate(i,j,tab_transition.GetValue(i,j));
		}
	}
	// Setup diffusion data
	int zpos = 0;
	assert(tab_layer_structure.NumRows() == NUM_LAYERS);
	for (int row = 0; row < NUM_LAYERS; row++)
	{
		double D = tab_layer_structure.GetValue(row,0);
		for (int layer = 0; layer < (int)(LayerThicknessFraction[row]*(double)(nk))+1 && zpos < nk; layer++)
		{
			for (int i = 0; i < ni; i++)
			{
				for (int j = 0; j < nj; j++)
				{
					des_model->setDiffusionConst(D,i,j,zpos);
				}
			}
			zpos++;
		}
	}
	assert(zpos == nk);
	// Initialize cells
	int xstart = gsl_rng_get(r)%ni;
	int ystart = gsl_rng_get(r)%nj;
	int zstart = gsl_rng_get(r)%nk;
	for(int i = 0; i < ni; i++)
	{
		for(int j = 0; j < nj; j++)
		{
			for(int k = 0; k < nk; k++)
			{
				if (i == xstart && j == ystart && k == zstart)
					des_model->addInitialCell(i,j,k,2);
				else
					des_model->addInitialCell(i,j,k,1);
			}
		}
	}
	// Create the simulator
	sim = new des::Simulator(des_model);
	// Done!
	return true;
}

int main(int argc, char **argv)
{
	double tStop = 600.0, bOutputPeriod = tStop/10.0;

	for( int i = 1; i < argc; ++i )
    {
		if( strcmp( argv[i], "-var" ) == 0 && ++i < argc )
        {
            InputFileName = argv[i];
        }
		else if( strcmp( argv[i], "-output_period" ) == 0 && ++i < argc)
        {
			bOutputPeriod = atof(argv[i]);
		}
		else if( strcmp( argv[i], "-end_time" ) == 0 && ++i < argc )
        {
			tStop = atof( argv[i] );
        }
		else if( strcmp( argv[i], "-ranseed" ) == 0 && ++i < argc )
        {
            ranseed = (unsigned)atol( argv[i] );
        }
	}
	//Initialize the variables and tables 
	InitFunction();
	InitModelCondition();
	// Run the simulation
	int seq_num = 0;
	double tL = 0.0;
	while (sim->nextEventTime() < tStop)
	{
		if ((int)(sim->nextEventTime()/bOutputPeriod) >= seq_num)
		{	
			PrintCSV(seq_num++,tL);
		}
		tL = sim->nextEventTime();
		sim->execNextEvent();
	}
	PrintCSV(seq_num++,tL);
	// Cleanup
	PreKill();
	return 0;
}
