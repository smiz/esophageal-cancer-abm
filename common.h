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
		 * Select direction at random.
		 */
		void direction(int& dx, int& dy, int& dz);
		/**
		 * Get the size if mm of a cell side.
		 */
		double cell_size() const { return dx; }
		/**
		 * Set this before you do anything else.
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
		void set_mutations_per_year(double events_per_year, int type) {
			mutate_time[type] = stem_cells_per_mm2*dx*dx/events_per_year;
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
		double set_diffusion_rate(double mm2_year) {
			diff_time = (dx*dx)/mm2_year;
		}
		/**
		 * Wrap the point into the cellspace or return false if the
		 * point is not in the wrapped cell space.
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
		 * Load parameter data from a text file.
		 */
		void load_from_file(const char* filename);

		static void deleteInstance();
	private:
		Parameters();
		~Parameters();
		gsl_rng *r; // RNG and Distribution package
		int nx, ny, nz; // Number of cells in each direction
		double dx; // Size of a cell
		double diff_time; // Mean time to a diffusion event
		double mutate_time[NUM_CELL_TYPES];
		double stem_cells_per_mm2;
		static Parameters* inst; // The singleton instance
};

#endif