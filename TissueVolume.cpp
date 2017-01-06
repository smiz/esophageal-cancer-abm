#include "TissueVolume.h"

TissueVolume::TissueVolume(int iType, int x, int y, int z):
	adevs::Atomic<adevs::CellEvent<int> >(),
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
	return (ttm < tte) ? ttm : tte;
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
		// Otherwise never mutate
		else ttm = adevs_inf<double>();
		// If we can expand now then do so
		if (iType >= DYSPLASIA)
			tte = p->exponential(p->get_expand_interval());
		// Otherwise never expand
		else tte = adevs_inf<double>();
	}
	// Otherwise we will move
	else 
	{
		// Better be expanding if we didn't mutate
		assert(iType == DYSPLASIA || iType == CANCER);
		// Reduce our time to mutate
		if (ttm < adevs_inf<double>()) ttm -= tte;
		assert(iType == CANCER || ttm < adevs_inf<double>());
		// Expand again 
		tte = p->exponential(p->get_expand_interval());
	}
	assert(tte >= 0.0);
	assert(ttm >= 0.0);
}

void TissueVolume::delta_ext(double e, const adevs::Bag<adevs::CellEvent<int> >& xb)
{
	// Reduce the time to next event
	if (ttm < adevs_inf<double>()) ttm -= e;
	if (tte < adevs_inf<double>()) tte -= e;
	// Get the most aggressive type that is visit us
	int newType = 0;
	for (auto x : xb)
	{
		assert(x.value == DYSPLASIA || x.value == CANCER);
		if (x.value > newType) newType = x.value;
	}
	// If we are going to be invaded, reset the event timers
	if
	(
		newType != iType && // No need to change in this case
		(
			newType == CANCER // Cancer always spreads
				||
			(iType == BE && newType == DYSPLASIA) // Dysplasia can spread in BE
		)
	)
	{
		iType = newType;
		Parameters* p = Parameters::getInstance();
		// Expand
		tte = p->exponential(p->get_expand_interval());
		// Mutate if we can
		if (p->get_mutation_interval(iType) < adevs_inf<double>())
			ttm = p->exponential(p->get_mutation_interval(iType));
		else
			ttm = adevs_inf<double>();
	}
}

void TissueVolume::delta_conf(const adevs::Bag<adevs::CellEvent<int> >& xb)
{
	delta_int();
	delta_ext(0.0,xb);
	assert(tte >= 0.0);
	assert(ttm >= 0.0);
}

void TissueVolume::output_func(adevs::Bag<adevs::CellEvent<int> >& yb)
{
	// Output our type if we are expanding
	if (tte < ttm)
	{
		assert(iType == CANCER || iType == DYSPLASIA);
		int dx = 0, dy = 0, dz = 0;
		adevs::CellEvent<int> out;
		// Cancer can spread anywhere
		if (iType == CANCER)
			Parameters::getInstance()->direction(dx,dy,dz);
		// Dysplasia is stuck on the surface
		else
			Parameters::getInstance()->direction(dx,dy);
		dx += x; dy += y; dz += z;
		// If direction is out of the space, then no output
		if (!Parameters::getInstance()->wrap(dx,dy,dz))
			return;
		out.x = dx; out.y = dy; out.z = dz;
		out.value = iType;
		yb.insert(out);
	}
}

