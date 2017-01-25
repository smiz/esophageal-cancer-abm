#ifndef _common_h_
#define _common_h_
#include "adevs.h"
#include <gsl/gsl_randist.h>

/**
 * Types of cells in the model. These
 * must be in order of evolution.
 */
#define NORMAL 0
#define BE 1
#define DYSPLASIA 2
#define CANCER 3
#define NUM_CELL_TYPES 4

/**
 * This class stores all globally visible model parameters.
 */
class Parameters
{
	public:
		/**
		 * Set the random number seed.
		 */
		void set_seed(unsigned long seed);
		/**
		 * Use the global random number generator to
		 * sample an exponential distribution with mean mu.
		 */
		double exponential(double mu);
		/**
		 * Return a number in [0,1]
		 */
		double uniform();
		/** Sample a normal distribution */
		double normal(double mean, double std_dev);
		/**
		 * Select a 3D direction at random.
		 */
		void direction(int& dx, int& dy, int& dz);
		/**
		 * Select a 2D direction at random.
		 */
		void direction(int& dx, int& dy);
		/**
		 * Get the size in mm of a tissue volume block.
		 */
		double cell_size() const { return dx; }
		/**
		 * Set the tissue volume size before you do anything else.
		 */
		void cell_size(double ddxx) { dx = ddxx; }
		/**
		 * Get or set the number of cells in the x, y, and z
		 * directions.
		 */
		int xdim() const { return nx; }
		int ydim() const { return ny; }
		int zdim() const { return nz; }
		void xdim(int x) { nx = x; }
		void ydim(int y) { ny = y; }
		void zdim(int z) { nz = z; }
		/**
		 * Set the mutation interval for a single stem cell. Make sure
		 * dx is set properly first.
		 */
		void set_mutations_per_year(double mutations, int type) {
			mutate_time[type] = 1.0/(stem_cells_per_mm2*dx*dx*mutations);
		}
		/**
		 * Get the mutation interval (in years) for a cell site in the model
		 */
		double get_mutation_interval(int type) const { return mutate_time[type]; }
		/**
		 * Get the time (in years) until a dysplasia or cancer cell will expand.
		 */
		double get_expand_interval() const { return diff_time; }
		/**
		 * Set the diffusion rate for malignant cells.
		 */
		void set_diffusion_rate(double mm2_year) {
			diff_time = (dx*dx)/mm2_year;
		}
		/**
		 * Wrap the x,y,z point into the cellspace or return false if the
		 * point is not in the wrapped space.
		 */
		bool wrap(int& x, int& y, int& z) const;
		/**
		 * Get the singleton instance of this object.
		 */
		static Parameters* getInstance();
		/**
		 * Set stem cells per square mm. Set this after dx.
		 */
		void set_stem_cells_per_mm2(double count) {
			stem_cells_per_mm2 = count;
		}
		/**
		 * Set the mean age for onset of BE
		 */
		void be_onset_age(double years) {
			be_onset = years;
		}
		/**
		 * Get the mean age for onset of BE
		 */
		double be_onset_age() const {
			return be_onset;
		}
		/**
		 * Load parameter data from a text file.
		 */
		void load_from_file(const char* filename);
		/// Delete the current singleton instance
		static void deleteInstance();
	private:
		/// This is a singleton, so constructor and destructor are private
		Parameters();
		Parameters(const Parameters&){}
		Parameters& operator=(const Parameters& other) { return *this; } 
		~Parameters();
		gsl_rng *r; // RNG and Distribution package
		int nx, ny, nz; // Number of cells in each direction
		double dx; // Size of a cell
		double diff_time; // Mean time to a diffusion event
		double mutate_time[NUM_CELL_TYPES];
		double stem_cells_per_mm2;
		double be_onset;
		static Parameters* inst; // The singleton instance
};

#endif
