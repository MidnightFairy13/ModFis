#pragma once
#include "Null_Cond.h"
#include <iostream>
/*
При записи элементов стоит учесть:
последние элементы MET и GAS будут расчитыватся относительно предыдущих

Material (name, reag, Nu, Y.Bound);
Значения Nu ставятся в соответсвии с уравнением
#пример
2 Al + 3 CO2 = 1 Al2O3 + 3 CO + 0 N2
reag:
1: реагирующее вещество
-1: продукт реакции
0: нереагирующее вещество

для газа
U = { Ro, RoU, RoV, RoE, Ro_i... }
F = { RoUn, RoUn^2+P, RoUnUt, RoUnH, Ro_iUn... }

для частиц
U = { Ns, Ro, RoU, RoV, RoTCs, Ro_i... }
F = { NsUn, RoUn, RoUn^2, RoUnUt, RoUnTCs, Ro_iUn... }

pos_flow_ - номера потоков для вывода на консоль
*/
Inf_Cond::Inf_Cond()
{

	Gas Ox("CO2", 1, 3, 0.3), Pr_g("CO", -1, 3, 1 - Ox.Bound_Y), Inert("N2", 0);
	GS = { Ox, Pr_g, Inert };
	Meas_Gas();
	eq_fi = (int)GS.size() - 1;
	pos_flows_fluid = { 0, 3 };
	if (!part) in.z = 0;
	if (part)
	{
		Metal R("Al", 1, 2, 1), Pr_m("Al2O3", -1, 1, 1 - R.Bound_Y);
		MT = { R, Pr_m };
		Meas_Metal();
		eq_p = 5, eq_pi = (int)MT.size() - 1;
		pos_flows_part = { 1, 4 };
		if (burn) pos_flows_part.push_back(5);
	}
	if (eq_fi == -1 || eq_pi == -1) std::cerr << "Error eq_i < 0", exit(1);
}
void Inf_Cond::Meas_Metal()
{
	Meas meas;
	double Sum = 0;
	for (int pos_mat = 0; pos_mat < MT.size(); pos_mat++)
	{
		MT[pos_mat].Cs /= meas.R;
		MT[pos_mat].Ro /= meas.Ro;
		Sum += MT[pos_mat].Bound_Y / MT[pos_mat].Ro;
	}
	in.M1 = meas.pi * pow(in.D, 3) / 6. / Sum;
}
void Inf_Cond::Meas_Gas()
{
	Meas meas;
	for (int pos_mat = 0; pos_mat < GS.size(); pos_mat++)
	{
		GS[pos_mat].Cp /= meas.R, GS[pos_mat].Cv /= meas.R,
			GS[pos_mat].Lam /= meas.Lam, GS[pos_mat].Mu /= meas.Mu, GS[pos_mat].R /= meas.R;
	}
}

Meas::Meas()
{
	D = 1e+6 * X;
	V = pow(P / Ro, 0.5), Time = X / V;
	T = V * V / R;
	Mu = P * Time, Lam = R * V * X * Ro;
	G = Ro * V, Q = V * V;
}

Null_Fluid::Null_Fluid()
{
	Meas meas;
	U /= meas.V, V /= meas.V, P /= meas.P, Ro /= meas.Ro;
}

Null_Part::Null_Part()
{
	Meas meas;
	U /= meas.V, V /= meas.V, T /= meas.T, Ro /= meas.Ro;
}

Inlet::Inlet()
{
	Meas meas;
	G /= -meas.G, T /= meas.T;
	P /= meas.P, Ro / meas.Ro;
	D /= meas.D, Q /= meas.Q;
	Rot /= meas.Ro, U01 /= -meas.V, P1 /= pow(meas.P, nu);
}

Outlet::Outlet()
{
	Meas meas;
	P /= meas.P, Ro /= meas.Ro;
}
