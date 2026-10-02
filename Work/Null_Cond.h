#pragma once
#include <iomanip>
#include <vector>
#include "Material.h"
struct Meas
{
	double X = 1., R = 287., Ro = 1., P = 1e+5, pi = 3.1415926535, eps = 1e-4,
		D, V, Time, T, Mu, Lam, G, Q;
	Meas();
};
struct Inlet
{
	double G = 100., T = 1300.,
		Rot = 1700, nu = 0.4, U01 = 5e-3, P1 = pow(101325, nu),
		P = 1e+5, M = 0., Ro = 1.,
		D = 20., z = 0.1, Q = 20e+6, M1 = 0.;
	Inlet();
};
struct Outlet
{
	double P = 1e+5, Ro = 1.;
	Outlet();
};
struct Null_Fluid
{
	double U = 0., V = 0., P = 1e+5, Ro = 1.;
	Null_Fluid();
};
struct Null_Part
{
	double U = 0., V = 0., T = 0., Ro = 0., Ns = 0.;
	Null_Part();
};
/*
Metod_Mdt
1:Беляев-Фролов
2:Бекстед

Metod_axis
1:С правыми частями
2:Без правых частей

Metod_mass_flow
1: G, T
2: Rot, U01, Nu, Tr

set:
0: до установления
n: количество циклов расчета

ddi - количество расчетов в одном цикле
diter - через сколько расчетов вывести на консоль потоки

*/
class Inf_Cond
{
public:
	int eq_f = 4, eq_p = 0, eq_fi = 0, eq_pi = 0,
		Metod_Mdt = 1, Metod_axis = 2, Metod_mass_flow = 2, set = 0, ddi = 100, diter = ddi * ddi;
	double ddt = 0.8;
	bool part = 1, burn = 1, continue_file = 0;
	Inlet in;
	Outlet out;
	std::vector<int> pos_flows_fluid, pos_flows_part;
	std::vector<Metal> MT;
	std::vector<Gas> GS;
	Inf_Cond();
private:
	void Meas_Metal();
	void Meas_Gas();
};