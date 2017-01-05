#ifndef _tissue_h_
#define _tissue_h_
#include "common.h"
#include "TissueVolume.h"

/**
 * Model of the esophagaes
 */
class Tissue:
	public adevs::CellSpace<int>::Cell
{
	public:
		TissueVolume(int iType, int x, int y, int z);
		double ta();
		void delta_int();
		void delta_ext(double e, const adevs::Bag<adevs::CellEvent<int> >& xb);
		void delta_conf(const adevs::Bag<adevs::CellEvent<int> >& xb);
		void output_func(adevs::Bag<adevs::CellEvent<int> >& yb);
		int xpos() const { return x; }
		int ypos() const { return y; }
		int zpos() const { return z; }
		int itype() const { return iType; }
		~TissueVolume();
	private:
		int iType; // Type of cell
		double ttm, tte; // Time to mutate and expand
		const int x, y, z;
};

#endif

