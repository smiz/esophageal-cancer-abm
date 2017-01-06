#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "common.h"
#include "TissueVolume.h"
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
// Fractional thickness of each layer
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

// The TissueVolume objects in this CellSpace comprise the dynamic part of the model
static CellSpace<int>* tissue;
// The Simulator handles time management, etc.
static Simulator<CellEvent<int> >* sim;
// File for parameters not hard coded into the model.
static std::string inputData = "input.txt";
// Time points for taking biopsies
static list<double> biopsy;
// Onset age for be
static double be_onset;

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
	int types[NUM_CELL_TYPES] = { 0, 0, 0, 0 };
	char filename[100];
	sprintf(filename,"tumor.csv.%d",seq_num);
	// Output the paraview file format
	ofstream fout(filename);
	fout << "xcoord,ycoord,zcoord,type" << endl;
	for (int i = 0; i < ni; i++)
		for (int j = 0; j < nj; j++)
			for (int k = 0; k < nk; k++)
			{
				int type =
					dynamic_cast<TissueVolume*>(tissue->getModel(i,j,k))->itype();
				types[type]++;
				if (type > BE)
				{
					fout << i << "," << j << "," << k << "," << type << endl;
				}
			}
	fout.close();
	// Report the time and counts of each cell type
	cout << "t = " << t << endl;
	for (int i = 0; i < NUM_CELL_TYPES; i++)
	{
		cout << names[i] << " : " << types[i] << " " << endl;
	}
}

//===========================================================================//
void InitModel(void)
{
	// Set the size of the simulation grid
	Parameters::getInstance()->cell_size(grid_size);
	Parameters::getInstance()->xdim(ni);
	Parameters::getInstance()->ydim(nj);
	Parameters::getInstance()->zdim(nk);
	// Load the free parameters
	Parameters::getInstance()->load_from_file(inputData.c_str());
	// Create the simulation grid
	tissue = new CellSpace<int>(ni,nj,nk);
	// Populate it with TissueVolume objects
	for (int i = 0; i < ni; i++)
	{
		for (int j = 0; j < nj; j++)
		{
			for (int k = 0; k < nk; k++)
			{
				// Top layer has BE
				if (k == 0)
					tissue->add(new TissueVolume(BE,i,j,k),i,j,k);
				// Everything else is initially normal
				else
					tissue->add(new TissueVolume(NORMAL,i,j,k),i,j,k);
			}
		}
	}
	// Get the BE onset age
	double mean_onset = Parameters::getInstance()->be_onset_age();
	be_onset = Parameters::getInstance()->exponential(mean_onset);
	// Sort the biopsies by age
	biopsy.sort();
	// Create the simulator for our tissue model
	sim = new Simulator<CellEvent<int> >(tissue);
}

int main(int argc, char **argv)
{
	// Read and apply the command line arguments
	for( int i = 1; i < argc; ++i )
	{
		if(strcmp( argv[i],"-var") == 0 && ++i < argc)
		{
			inputData = argv[i];
		}
		else if( strcmp( argv[i], "-ranseed" ) == 0 && ++i < argc )
		{
			unsigned ranseed = (unsigned)atol( argv[i] );
			Parameters::getInstance()->set_seed(ranseed);
		}
		else  
		{
			errno = 0;
			double age = strtod(argv[i],NULL);
			if (age <= 0.0 || errno != 0)
			{
				cout << "Illegal biopsy age " << argv[i] << endl;
				return 0;
			}
			biopsy.push_back(age);
		}
	}
	// Setup the model
	InitModel();
	// Run the simulation
	int seq_num = 0;
	while (!biopsy.empty())
	{
		if (sim->nextEventTime()+be_onset > biopsy.front())
		{
			PrintCSV(seq_num++,biopsy.front());
			biopsy.pop_front();
		}
		else
			sim->execNextEvent();
	}
	// Cleanup
	delete sim;
	delete tissue;
	Parameters::deleteInstance();
	return 0;
}
