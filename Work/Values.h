#pragma once
#include "Fisics_El.h"
#include "Null_Cond.h"
namespace val
{
	void plus(std::vector<double>& f, std::vector<double> df);
	void mult(std::vector<double>& f, double a);
	double Nu_M(Inf_Cond& inf, int pos_mat, bool fluid = 1);
	double Ros_0(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Ds(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Dsi(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat);
	double N_gas(Inf_Cond& inf, gm::OPT_CELL& cell);
	double N_i(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat);
	double C_mat(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat);
	double X_eff(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Cp(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Cp(Inf_Cond& inf);
	double Cv(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Cv(Inf_Cond& inf);
	double R(Inf_Cond& inf, gm::OPT_CELL& cell);
	double R(Inf_Cond& inf);
	double T(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Cs(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Cs(Inf_Cond& inf);
	double Mu(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Lam(Inf_Cond& inf, gm::OPT_CELL& cell);
	double k(Inf_Cond& inf, gm::OPT_CELL& cell);
	double k(Inf_Cond& inf);
	double H(Inf_Cond& inf, gm::OPT_CELL& cell);
	double M(Inf_Cond& inf, gm::OPT_CELL& cell, double U);
	double M(Inf_Cond& inf, gm::OPT_CELL& cell);
	double C(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Re(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Pr(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Nu(Inf_Cond& inf, gm::OPT_CELL& cell);
	double fd(Inf_Cond& inf, gm::OPT_CELL& cell);
	double fi_1(Inf_Cond& inf, gm::OPT_CELL& cell);
	double fi_2(Inf_Cond& inf, gm::OPT_CELL& cell);
	double Mdt(Inf_Cond& inf, gm::OPT_CELL& cell);
}
