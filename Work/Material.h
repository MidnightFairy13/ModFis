#pragma once
#include <string>
#include <vector>

struct Metal
{
	std::string name;
	int reag;
	double Nu, Bound_Y, Cs, Ro, M;
	Metal();
	Metal(std::string _name, int _reag, double _Nu = 0, double _Bound_Y = 0);
};
struct Gas
{
	std::string name;
	int reag;
	double Nu, Bound_Y, Cp, Cv, M, Mu, Lam, R;
	Gas();
	Gas(std::string _name, int _reag, double _Nu = 0, double _Bound_Y = 0);
};
struct Real_Metals
{
	std::vector<std::string> name;
	std::vector<double> Cs, Ro, M;
	Real_Metals();
	int POS_METAL(std::string _name);
};
struct Real_Gases
{
	std::vector<std::string> name;
	std::vector<double> Cp, Cv, M, Mu, Lam, R;
	Real_Gases();
	int POS_GAS(std::string _name);
};

