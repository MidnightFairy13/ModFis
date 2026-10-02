#pragma once
#include "Values.h"
#include <iostream>
namespace val
{
	void plus(std::vector<double>& f, std::vector<double> df)
	{
		for (int i = 0; i < f.size(); i++) f[i] += df[i];
	}

	void mult(std::vector<double>& f, double a)
	{
		for (int i = 0; i < f.size(); i++) f[i] *= a;
	}

	double Nu_M(Inf_Cond& inf, int pos_mat, bool fluid)
	{
		if (fluid) return inf.GS[pos_mat].reag * inf.GS[pos_mat].Nu * inf.GS[pos_mat].M;
		else return inf.MT[pos_mat].reag * inf.MT[pos_mat].Nu * inf.MT[pos_mat].M;
	}
	double Ros_0(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Ro = 0;
		for (int i = 0; i < inf.MT.size(); i++) Ro += cell.Part.Opt.Y[i] / inf.MT[i].Ro;
		return Ro;
	}
	double Ds(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		Meas meas;
		double D = 0;
		if (cell.Part.Opt.Ns != 0)
		{
			double T = cell.Part.Opt.Ro / cell.Part.Opt.Ns;
			if (Ros_0(inf, cell) < 0) return 0;
			D = pow(6. / meas.pi * T * Ros_0(inf, cell), 1. / 3.);
		}
		return D;
	}
	double Dsi(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat)
	{
		Meas meas;
		double D = 0;
		if (cell.Part.Opt.Ns != 0)
		{
			double T = cell.Part.Opt.Ro / cell.Part.Opt.Ns;
			if (cell.Part.Opt.Y[pos_mat] < 0) return 0;
			D = pow(6. / meas.pi * T * cell.Part.Opt.Y[pos_mat] / inf.MT[pos_mat].Ro, 1. / 3.);
		}
		return D;
	}
	double N_gas(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Ro * cell.Fluid.Opt.Y[i] / inf.GS[i].M;
		return Sum;
	}

	double N_i(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat)
	{
		return cell.Fluid.Opt.Ro * cell.Fluid.Opt.Y[pos_mat] / inf.GS[pos_mat].M;
	}
	double C_mat(Inf_Cond& inf, gm::OPT_CELL& cell, int pos_mat)
	{
		return N_i(inf, cell, pos_mat) / N_gas(inf, cell);
	}
	double X_eff(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		int O2 = -1, H2O = -1, CO2 = -1;
		double Xeff = 0;
		for (int i = 0; i <= inf.eq_fi; i++)
		{
			if (inf.GS[i].name == "O2") O2 = i;
			if (inf.GS[i].name == "H2O") H2O = i;
			if (inf.GS[i].name == "CO2") CO2 = i;
		}
		if (O2 != -1) Xeff += C_mat(inf, cell, O2);
		if (H2O != -1) Xeff += 0.6 * C_mat(inf, cell, H2O);
		if (CO2 != -1) Xeff += 0.22 * C_mat(inf, cell, CO2);
		return Xeff;
	}
	double Cp(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Y[i] * inf.GS[i].Cp;
		return Sum;
	}

	double Cp(Inf_Cond& inf)
	{
		double Sum = 0;
		for (int i = 0; i <= inf.eq_fi; i++) Sum += inf.GS[i].Bound_Y * inf.GS[i].Cp;
		return Sum;
	}

	double Cv(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Y[i] * inf.GS[i].Cv;
		return Sum;
	}
	double Cv(Inf_Cond& inf)
	{
		double Sum = 0;
		for (int i = 0; i <= inf.eq_fi; i++) Sum += inf.GS[i].Bound_Y * inf.GS[i].Cv;
		return Sum;
	}
	double R(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Y[i] * inf.GS[i].R;
		return Sum;
	}
	double R(Inf_Cond& inf)
	{
		double Sum = 0;
		for (int i = 0; i <= inf.eq_fi; i++) Sum += inf.GS[i].Bound_Y * inf.GS[i].R;
		return Sum;
	}
	double T(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return cell.Fluid.Opt.P / (cell.Fluid.Opt.Ro * R(inf, cell));
	}
	double Cs(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Part.Opt.Y.size(); i++) Sum += cell.Part.Opt.Y[i] * inf.MT[i].Cs;
		return Sum;
	}

	double Cs(Inf_Cond& inf)
	{
		double Sum = 0;
		for (int i = 0; i <= inf.eq_pi; i++) Sum += inf.MT[i].Bound_Y * inf.MT[i].Cs;
		return Sum;
	}

	double Mu(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Y[i] * inf.GS[i].Mu;
		return Sum;
	}

	double Lam(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Sum = 0;
		for (int i = 0; i < cell.Fluid.Opt.Y.size(); i++) Sum += cell.Fluid.Opt.Y[i] * inf.GS[i].Lam;
		return Sum;
	}
	double k(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return Cp(inf, cell) / Cv(inf, cell);
	}
	double k(Inf_Cond& inf)
	{
		return Cp(inf) / Cv(inf);
	}
	double H(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return k(inf, cell) / (k(inf, cell) - 1.) * cell.Fluid.Opt.P / cell.Fluid.Opt.Ro + (cell.Fluid.Opt.U * cell.Fluid.Opt.U + cell.Fluid.Opt.V * cell.Fluid.Opt.V) / 2.;
	}
	double M(Inf_Cond& inf, gm::OPT_CELL& cell, double U)
	{
		return U / C(inf, cell);
	}
	double M(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return cell.Fluid.Opt.U / C(inf, cell);
	}
	double C(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return sqrt(abs(k(inf, cell) * cell.Fluid.Opt.P / cell.Fluid.Opt.Ro));
	}
	double Re(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double U, V, d = 0;
		U = cell.Part.Opt.U - cell.Fluid.Opt.U;
		V = cell.Part.Opt.V - cell.Fluid.Opt.V;
		if (inf.part) d = Ds(inf, cell);
		return cell.Fluid.Opt.Ro * d * pow(U * U + V * V, 1 / 2) / Mu(inf, cell);
	}

	double Pr(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return Cp(inf, cell) * Mu(inf, cell) / Lam(inf, cell);
	}

	double Nu(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return 2. + 0.6 * pow(Re(inf, cell), 1 / 2.) * pow(Pr(inf, cell), 1 / 3.);
	}

	double fd(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		return 1. + 0.15 * pow(Re(inf, cell), 0.687);
	}

	double fi_1(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		double Fi;
		if (Ds(inf, cell) == 0 || Ros_0(inf, cell) < 0) return 0;
		else Fi = 18. * fd(inf, cell) * Mu(inf, cell) / (Ds(inf, cell) * Ds(inf, cell)) * Ros_0(inf, cell);
		if (Fi != Fi)
		{
			std::cout << "fd = " << fd(inf, cell) << std::endl
				<< "Mu = " << Mu(inf, cell) << std::endl
				<< "Ds = " << Ds(inf, cell) << std::endl
				<< "Ros_0 = " << Ros_0(inf, cell) << std::endl;
		}
		return Fi;
	}
	double fi_2(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		if (Cs(inf, cell) == 0) return 0;
		else return fi_1(inf, cell) * Nu(inf, cell) * Cp(inf, cell) / (3. * Pr(inf, cell) * Cs(inf, cell) * fd(inf, cell));
	}

	double Mdt(Inf_Cond& inf, gm::OPT_CELL& cell)
	{
		Meas meas;
		switch (inf.Metod_Mdt)
		{
		case 1:
		{
			double Tet = C_mat(inf, cell, 0), mdt, K;
			K = pow(10, 4.8) * meas.Time / (3. * 0.67 * pow(meas.D, 1.5));
			mdt = -meas.pi * inf.MT[0].Ro * pow(Dsi(inf, cell, 0), 1.5) * pow(Tet, 0.9) * K;
			return mdt;
		}
		case 2:
		{
			double Xeff = X_eff(inf, cell), mdt, K, a = 0.0244, n = 1.5;
			K = pow(10, 2.5) * pow(meas.P, 0.1) * pow(meas.T, 0.2) * meas.Time / (a * 2. * n * pow(meas.D, n));
			mdt = -meas.pi * inf.MT[0].Ro * pow(Dsi(inf, cell, 0), 3 - n) * Xeff * pow(cell.Fluid.Opt.P, 0.1) * pow(T(inf, cell), 0.2) * K;
			return mdt;
		}
		default: return 0;
		}
	}
}
