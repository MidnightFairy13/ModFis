#pragma once
#include <vector>
namespace gm
{
	struct FLUID
	{
		double Ro = 0., P = 0., U = 0., V = 0.;
		std::vector <double> Y;
	};
	struct PARTICLE
	{
		double  Ro = 0., T = 0., U = 0., V = 0., Ns = 0.;
		std::vector <double> Y;
	};
	struct OPT_FLUID
	{
		FLUID Opt;
		std::vector <double> Flow;
	};
	struct OPT_PART
	{
		PARTICLE Opt;
		std::vector <double> Flow;
	};
	struct OPT_CELL
	{
		OPT_FLUID Fluid;
		OPT_PART Part;
	};
	struct OPT_FACE
	{
		std::vector <double> Fluid_Flow, Part_Flow;
	};
}