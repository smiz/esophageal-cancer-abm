#ifndef _cell_h_
#define _cell_h_
#include "common.h"

/**
 * Model of a volume of tissue
 */
class TissueVolume:
	public adevs::Atomic<adevs::CellEvent<int> >
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
		void gc_output(adevs::Bag<adevs::CellEvent<int> >&){}
		~TissueVolume(){}
	private:
		int iType; // Type of cell
		double ttm, tte; // Time to mutate and expand
		const int x, y, z;
};

#endif

