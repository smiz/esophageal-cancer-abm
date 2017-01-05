#include "TissueVolume.h"

TissueVolume::TissueVolume(int iType, int x, int y, int z):
	adevs::CellSpace<int>::Cell(),
	iType(iType),
	ttm(adevs_inf<double>()),
	tte(adevs_inf<double>()),
	x(x),y(y),z(z)
{
	Parameters* p = Parameters::getInstance();
	// Only dysplasia and cancer can expand
	if (iType == DYSPLASIA || iType == CANCER)
		tte = p->exponential(p->get_expand_interval());
	// Anything might mutate
	if (p->get_mutation_interval(iType) < adevs_inf<double>())
		ttm = p->exponential(p->get_mutation_interval(iType));
}

double TissueVolume::ta()
{
	(ttm < tte) ? ttm : tte;
}

void TissueVolume::delta_int()
{
	Parameters* p = Parameters::getInstance();
	// We will mutate
	if (ttm < tte)
	{
		// Reduce our time to expand by what has passed
		if (tte < adevs_inf<double>()) tte -= ttm;
		// Cancer never mutates
		assert(iType < CANCER);
		// Advance our type
		iType++;
		// If we can mutate as the new type, pick a time to mutate
		if (p->get_mutation_interval(iType) < adevs_inf<double>())
			ttm = p->exponential(p->get_mutation_interval(iType));
		// If we can expand now then do so
		if (iType == DYSPLASIA)
		{
			assert(tte == adevs_inf<double>());
			tte = p->exponential(p->get_expand_interval());
		}
	}
	// Otherwise we will move
	else 
	{
		// Better be expanding if we didn't mutate
		assert(iType == DYSPLASIA || iType == CANCER);
		// Reduce our time to mutate
		if (ttm < adevs_inf<double>()) ttm -= tte;
		// Expand again if we can
		tte = p->exponential(p->get_expand_interval());
	}
}

void TissueVolume::delta_ext(double e, const adevs::Bag<adevs::CellEvent<int> >& xb)
{
	// Reduce the time to next event
	if (ttm < adevs_inf<double>()) ttm -= e;
	if (tte < adevs_inf<double>()) tte -= e;
	// Get the most aggressive type that is visit us
	int newType = 0;
	for (auto x : xb)
		if (x.value > newType) newType = x.value;
	// If we are going to be invaded, reset the event timers
	if (newType > iType)
	{
		iType = newType;
		Parameters* p = Parameters::getInstance();
		if (iType == DYSPLASIA || iType == CANCER)
			tte = p->exponential(p->get_expand_interval());
		if (p->get_mutation_interval(iType) < adevs_inf<double>())
			ttm = p->exponential(p->get_mutation_interval(iType));
	}
}

void TissueVolume::delta_conf(const adevs::Bag<adevs::CellEvent<int> >& xb)
{
	delta_int();
	delta_ext(0.0,xb);
}

void TissueVolume::output_func(adevs::Bag<adevs::CellEvent<int> >& yb)
{
	// Output our type if we are expanding
	if (tte < ttm)
	{
		int dx, dy, dz;
		adevs::CellEvent<int> out;
		Parameters::getInstance()->direction(dx,dy,dz);
		dx += x; dy += y; dz += z;
		Parameters::getInstance()->wrap(dx,dy,dz);
		out.x = dx; out.y = dy; out.z = dz;
		out.value = iType;
		yb.insert(out);
	}
}

TissueVolume::~TissueVolume()
{
}

