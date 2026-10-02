#pragma once
#include "Material.h"
#include <iostream>
Real_Metals::Real_Metals()
{
	name = { "Al", "Al2O3" };
	Cs = { 920, 780 };
	Ro = { 2.7e+3, 3.99e+3 };
	M = { 27e-3, 102e-3 };
}
int Real_Metals::POS_METAL(std::string _name)
{
	int i;
	for (i = 0; i < name.size() && _name != name[i]; i++);
	if (i == name.size()) std::cerr << "Error: not found metal", exit(1);
	return i;
}
Real_Gases::Real_Gases()
{
	name = { "CO2", "CO", "N2" };
	Cp = { 837, 1040, 1050 };
	Cv = { 648, 743, 753 };
	Mu = { 149e-7, 166e-7, 165e-7 };
	Lam = { 166e-4, 232.6e-4, 0.026 };
	M = { 44e-3, 28e-3, 28e-3 };
	R = { 189, 297, 297 };
}
int Real_Gases::POS_GAS(std::string _name)
{
	int i;
	for (i = 0; i < name.size() && _name != name[i]; i++);
	if (i == name.size()) std::cerr << "Error: not found gas", exit(1);
	return i;
}
Metal::Metal()
{
	Bound_Y = 0, Nu = 0, reag = 0;
}
Metal::Metal(std::string _name, int _reag, double _Nu, double _Bound_Y)
{
	Real_Metals Met;
	int i = Met.POS_METAL(_name);
	name = _name;
	Cs = Met.Cs[i], Ro = Met.Ro[i], M = Met.M[i];
	Bound_Y = _Bound_Y, Nu = _Nu, reag = _reag;
}
Gas::Gas()
{
	Bound_Y = 0, Nu = 0, reag = 0;
}
Gas::Gas(std::string _name, int _reag, double _Nu, double _Bound_Y)
{
	Real_Gases Gas;
	int i = Gas.POS_GAS(_name);
	name = _name;
	Cp = Gas.Cp[i], R = Gas.R[i], M = Gas.M[i], Cv = Gas.Cv[i], Lam = Gas.Lam[i], Mu = Gas.Mu[i];
	Bound_Y = _Bound_Y, Nu = _Nu, reag = _reag;
}