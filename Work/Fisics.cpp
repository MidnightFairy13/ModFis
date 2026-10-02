#pragma once
#include <Windows.h>
#include <iostream>
#include "File.h"
#include <thread>
gm::FISICS::FISICS()
{
	std::cout << " Create Fisics Model...";
	Cells.resize(Mesh.count_all_cells);
	Faces.resize(Mesh.count_all_faces);
	for (int i = 0; i < Cells.size(); i++)
	{
		START_FLUID(Cells[i]);
		if (inf.part) START_PART(Cells[i]);
	}
	for (int i = 0; i < Faces.size(); i++)
	{
		Faces[i].Fluid_Flow.resize(inf.eq_f + inf.eq_fi);
		if (inf.part) Faces[i].Part_Flow.resize(inf.eq_p + inf.eq_pi);
	}
	std::cout << " OK" << std::endl;
}

void gm::FISICS::START_FLUID(gm::OPT_CELL& cell)
{
	cell.Fluid.Flow.resize(inf.eq_f + inf.eq_fi);
	Null_Fluid fl;
	int i = 0;
	cell.Fluid.Opt.Y.resize(inf.GS.size());
	for (; i < inf.eq_fi; i++) cell.Fluid.Opt.Y[i] = 0;
	cell.Fluid.Opt.Y[i] = 1;
	cell.Fluid.Opt.U = fl.U, cell.Fluid.Opt.V = fl.V, cell.Fluid.Opt.P = fl.P, cell.Fluid.Opt.Ro = fl.Ro;
	cell.Fluid.Flow[0] = cell.Fluid.Opt.Ro;
	cell.Fluid.Flow[1] = cell.Fluid.Opt.Ro * cell.Fluid.Opt.U;
	cell.Fluid.Flow[2] = cell.Fluid.Opt.Ro * cell.Fluid.Opt.V;
	cell.Fluid.Flow[3] = cell.Fluid.Opt.Ro * (cell.Fluid.Opt.U * cell.Fluid.Opt.U + cell.Fluid.Opt.V * cell.Fluid.Opt.V) / 2. + cell.Fluid.Opt.P / (val::k(inf, cell) - 1.);
	for (i = 0; i < inf.eq_fi; i++) cell.Fluid.Flow[inf.eq_f + i] = cell.Fluid.Flow[0] * cell.Fluid.Opt.Y[i];
}
void gm::FISICS::START_PART(gm::OPT_CELL& cell)
{
	cell.Part.Flow.resize(inf.eq_p + inf.eq_pi);
	Null_Part p;
	int i = 0;
	cell.Part.Opt.Y.resize(inf.MT.size());
	for (; i < cell.Part.Opt.Y.size(); i++) cell.Part.Opt.Y[i] = 0;
	cell.Part.Opt.U = p.U, cell.Part.Opt.V = p.V, cell.Part.Opt.Ro = p.Ro, cell.Part.Opt.T = p.T, cell.Part.Opt.Ns = p.Ns;
	cell.Part.Flow[0] = cell.Part.Opt.Ns;
	cell.Part.Flow[1] = cell.Part.Opt.Ro;
	cell.Part.Flow[2] = cell.Part.Opt.Ro * cell.Part.Opt.U;
	cell.Part.Flow[3] = cell.Part.Opt.Ro * cell.Part.Opt.V;
	cell.Part.Flow[4] = cell.Part.Opt.Ro * cell.Part.Opt.T * val::Cs(inf, cell);
	for (i = 0; i < inf.eq_pi; i++) cell.Part.Flow[inf.eq_p + i] = cell.Part.Flow[1] * cell.Part.Opt.Y[i];
}

void gm::FISICS::CALCULATING()
{
	std::string stop_string;
	std::vector<double> max, maxs;
	max.resize(inf.pos_flows_fluid.size()), maxs.resize(inf.pos_flows_part.size());
	bool stop = 0, fin = 0;
	double dt, time = 0;
	int finish = inf.set, count = 0;
	if (finish == 0) fin = 1;
	std::cout << " Calculation..." << std::endl;
	for (; !stop; stop_string.clear())
	{
		val::mult(max, 0), val::mult(maxs, 0);
		for (int i = 0; (finish > count || fin) && i < inf.ddi; i++)
		{
			dt = tau(), time += dt;
			FLOW_FACE();
			FLOW_CELL(dt);
			count++;
			if (count % inf.diter == 0)
			{
				std::cout << std::endl << "time: " << time << ", iter: " << count << std::endl;
				if (inf.burn) std::cout << std::endl << "burn: " << Tet() << std::endl;
			}
			stop = STOP();
			if (stop) break;
			stop_string += FLAG_FLOW(count, max);
			if (inf.part) stop_string += FLAG_FLOW(count, maxs, 0);
		}
		if (!stop) stop = WORK_STOP(count, finish, stop_string);
		if (finish != 0) fin = 0;
	}
}

std::string gm::FISICS::INF_FLOW(int pos_flow, bool fluid)
{
	switch (pos_flow)
	{
	case 0:
	{
		if (fluid) return "RoU";
		else return "NsU";
	}
	case 1:
	{
		if (fluid) return "RoU^2+P";
		else return "RosUs";
	}
	case 2:
	{
		if (fluid) return "RoVU";
		else return "RosUs^2";
	}
	case 3:
	{
		if (fluid) return "RoUH";
		else return "RosVsUs";
	}
	case 4:
	{
		if (fluid) return  "RoU(" + inf.GS[pos_flow - inf.eq_f].name + ')';
		else return "RosUsTsCs";
	}
	default:
	{
		if (fluid) return "RoU(" + inf.GS[pos_flow - inf.eq_f].name + ')';
		else return "RosUs(" + inf.MT[pos_flow - inf.eq_p].name + ')';
	}
	}
}

inline void gm::FISICS::MAX_DIF(int pos_cell, double& U, double& V)
{
	double u, v, c;
	c = val::C(inf, Cells[pos_cell]);
	u = Cells[pos_cell].Fluid.Opt.U, v = Cells[pos_cell].Fluid.Opt.V;
	U = abs(u + c), V = abs(v + c);
	if (U < abs(u - c)) U = abs(u - c);
	if (V < abs(v - c)) V = abs(v - c);
}

double gm::FISICS::tau()
{
	double tau = 100;
	double U = 0, V = 0, t_x = 0, t_y = 0;
	for (int i = 0; i < Mesh.count_all_cells; i++)
	{
		MAX_DIF(i, U, V);
		t_x = Mesh.cells[i].h / U;
		t_y = Mesh.cells[i].h / V;
		if (tau > t_x * t_y / (t_x + t_y)) tau = t_x * t_y / (t_x + t_y);
	}
	return tau * inf.ddt;
}

inline double gm::FISICS::Un(int pos_cell, int pos_face, bool fluid)
{
	double u = 0, v = 0;
	if (fluid) u = Cells[pos_cell].Fluid.Opt.U, v = Cells[pos_cell].Fluid.Opt.V;
	else u = Cells[pos_cell].Part.Opt.U, v = Cells[pos_cell].Part.Opt.V;
	return u * Mesh.faces[pos_face].n.x + v * Mesh.faces[pos_face].n.y;
}

inline double gm::FISICS::Ut(int pos_cell, int pos_face, bool fluid)
{
	double u = 0, v = 0;
	if (fluid) u = Cells[pos_cell].Fluid.Opt.U, v = Cells[pos_cell].Fluid.Opt.V;
	else u = Cells[pos_cell].Part.Opt.U, v = Cells[pos_cell].Part.Opt.V;
	return u * Mesh.faces[pos_face].t.x + v * Mesh.faces[pos_face].t.y;
}

inline void gm::FISICS::T_1(double& f1, double& f2, int pos_face)
{
	double _f1 = f1, _f2 = f2;
	f1 = _f1 * Mesh.faces[pos_face].n.x + _f2 * Mesh.faces[pos_face].t.x;
	f2 = _f1 * Mesh.faces[pos_face].n.y + _f2 * Mesh.faces[pos_face].t.y;
}

inline void gm::FISICS::T_1(int pos_face, bool fluid)
{
	if (fluid) T_1(Faces[pos_face].Fluid_Flow[1], Faces[pos_face].Fluid_Flow[2], pos_face);
	else T_1(Faces[pos_face].Part_Flow[2], Faces[pos_face].Part_Flow[3], pos_face);
}

void gm::FISICS::F_FLUID(int pos_face, double _U, double _V, double _P, double _Ro, bool inlet)
{
	int pos_cell = Mesh.faces[pos_face].cell_l;
	double k, Y;
	if (inlet) k = val::k(inf);
	else k = val::k(inf, Cells[pos_cell]);
	Faces[pos_face].Fluid_Flow[0] = _Ro * _U;
	Faces[pos_face].Fluid_Flow[1] = _Ro * _U * _U + _P;
	Faces[pos_face].Fluid_Flow[2] = _Ro * _U * _V;
	Faces[pos_face].Fluid_Flow[3] = _Ro * _U * ((_U * _U + _V * _V) / 2. + _P / _Ro * k / (k - 1.));
	T_1(pos_face);
	for (int j = 0; j < inf.eq_fi; j++)
	{
		if (inlet) Y = inf.GS[j].Bound_Y;
		else Y = Cells[pos_cell].Fluid.Opt.Y[j];
		Faces[pos_face].Fluid_Flow[inf.eq_f + j] = Faces[pos_face].Fluid_Flow[0] * Y;
	}
}

inline std::vector<double> gm::FISICS::F_FLUID(int pos_cell, int pos_face)
{
	std::vector<double> Fl(inf.eq_f + inf.eq_fi);
	double Ro, U, V, P;
	Ro = Cells[pos_cell].Fluid.Opt.Ro;
	P = Cells[pos_cell].Fluid.Opt.P;
	U = Un(pos_cell, pos_face);
	V = Ut(pos_cell, pos_face);
	Fl[0] = Ro * U;
	Fl[1] = Ro * U * U + P;
	Fl[2] = Ro * U * V;
	Fl[3] = Ro * U * ((U * U + V * V) / 2. + P / Ro * val::k(inf, Cells[pos_cell]) / (val::k(inf, Cells[pos_cell]) - 1.));
	for (int j = 0; j < inf.eq_fi; j++) Fl[inf.eq_f + j] = Fl[0] * Cells[pos_cell].Fluid.Opt.Y[j];
	return Fl;
}

inline void gm::FISICS::F_PART(int pos_face, double _U, double _V, double _T, double _Ro)
{
	int pos_cell = Mesh.faces[pos_face].cell_l;
	double Ros = _Ro * inf.in.z / (1. - inf.in.z), Ns = Ros / inf.in.M1;
	Faces[pos_face].Part_Flow[0] = Ns * _U;
	Faces[pos_face].Part_Flow[1] = Ros * _U;
	Faces[pos_face].Part_Flow[2] = Ros * _U * _U;
	Faces[pos_face].Part_Flow[3] = Ros * _U * _V;
	Faces[pos_face].Part_Flow[4] = Ros * _U * _T * val::Cs(inf);
	T_1(pos_face, 0);
	for (int j = 0; j < inf.eq_pi; j++) Faces[pos_face].Part_Flow[inf.eq_p + j] = Faces[pos_face].Part_Flow[1] * inf.MT[j].Bound_Y;
}

std::vector<double> gm::FISICS::F_PART(int pos_cell, int pos_face)
{
	std::vector<double> Pt(inf.eq_p + inf.eq_pi);
	double Ro, U, V, T, Ns;
	Ro = Cells[pos_cell].Part.Opt.Ro;
	T = Cells[pos_cell].Part.Opt.T;
	U = Un(pos_cell, pos_face, 0);
	V = Ut(pos_cell, pos_face, 0);
	Ns = Cells[pos_cell].Part.Opt.Ns;
	Pt[0] = Ns * U;
	Pt[1] = Ro * U;
	Pt[2] = Ro * U * U;
	Pt[3] = Ro * U * V;
	Pt[4] = Ro * U * T * val::Cs(inf, Cells[pos_cell]);
	T_1(Pt[2], Pt[3], pos_face);
	for (int j = 0; j < inf.eq_pi; j++) Pt[inf.eq_p + j] = Pt[1] * Cells[pos_cell].Part.Opt.Y[j];
	return Pt;
}

std::string gm::FISICS::FLAG_FLOW(int iter, std::vector<double>& max, bool fluid)
{
	std::vector<double> G;
	Meas meas;
	std::string flag;
	int eq = inf.eq_p;
	if (fluid) eq = inf.eq_f;
	std::vector<int> pos_flows;
	if (fluid) pos_flows = inf.pos_flows_fluid;
	else pos_flows = inf.pos_flows_part;
	G = dG(iter, fluid);
	for (int i = 0; i < G.size(); i++)
	{
		if (G[i] > max[i]) max[i] = G[i];
		if (iter % inf.ddi == 0)
		{
			if (pos_flows[i] < eq)
				if (max[i] < meas.eps) flag += INF_FLOW(pos_flows[i], fluid) + ' ';
		}
	}
	return flag;
}

inline bool gm::FISICS::WORK_STOP(int count, int& finish, std::string stop_string)
{
	bool stop = 0, fl = 0, fin = 0;
	fl = 0/*FLOW_STOP(stop_string)*/, fin = FINISH(count, finish);
	if (fin + fl > 0) stop = 1, WRITE();
	return stop;
}

inline bool gm::FISICS::STOP()
{
	bool stop = 0;
	if (GetAsyncKeyState(VK_CONTROL) < 0 && GetAsyncKeyState(0x51) < 0) std::cin.ignore(), stop = 1;
	if (stop)
	{
		std::cout << std::endl << "Stop? (0 or 1)", std::cin >> stop;
		WRITE();
	}
	return stop;
}

inline bool gm::FISICS::FINISH(int count, int& finish)
{
	bool flag = 0;
	if (count == finish)
	{
		int cont;
		std::cout << std::endl << " finish count " << std::endl;
		std::cout << " Continue? (input count) ", std::cin >> cont;
		if (cont > 0) { finish = cont + count; }
		else flag = 1;
	}
	return flag;
}

inline bool gm::FISICS::FLOW_STOP(std::string flow_stop)
{
	bool flag = 0;
	if (flow_stop.size() > 0)
	{
		std::cout << std::endl << "(flow_fluid < exp): " << flow_stop << std::endl;
		std::cout << " Stop? (0 or 1) ", std::cin >> flag;
	}
	return flag;
}

inline void gm::FISICS::WRITE()
{
	File f;
	bool write = 0;
	std::cout << " Writing? (input 0 or 1) ", std::cin >> write;
	if (write) f.Opt(*this);
}

double gm::FISICS::dGi(int pos_zone, int pos_flow, bool fluid)
{
	double F = 0;
	int i = Mesh.zone_faces[pos_zone].face_first, j = Mesh.zone_faces[pos_zone].face_last;
	for (int k = i; k <= j; k++) F += Fi(k, pos_flow, fluid);
	return F;
}
std::vector<double> gm::FISICS::dGi(int pos_zone, bool fluid)
{
	int size;
	if (fluid) size = (int)inf.pos_flows_fluid.size();
	else size = (int)inf.pos_flows_part.size();
	int i = Mesh.zone_faces[pos_zone].face_first, j = Mesh.zone_faces[pos_zone].face_last;
	std::vector<double> F(size);
	for (int k = i; k <= j; k++) val::plus(F, Fi(k, fluid));
	return F;
}
double gm::FISICS::Fi(int pos_face, int pos_flow, bool fluid)
{
	double F, S;
	if (fluid) F = Faces[pos_face].Fluid_Flow[pos_flow];
	else F = Faces[pos_face].Part_Flow[pos_flow];
	if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) S = Mesh.faces[pos_face].S * Mesh.faces[pos_face].r;
	else S = Mesh.faces[pos_face].S;
	return (F * Mesh.faces[pos_face].n.x + F * Mesh.faces[pos_face].n.y) * S;
}
std::vector<double> gm::FISICS::Fi(int pos_face, bool fluid)
{
	std::vector<int> pos_flows;
	if (fluid) pos_flows = inf.pos_flows_fluid;
	else pos_flows = inf.pos_flows_part;
	int size = (int)pos_flows.size();
	double S;
	std::vector<double> F(size);
	for (int pos = 0; pos < size; pos++)
	{
		if (fluid) F[pos] = Faces[pos_face].Fluid_Flow[pos_flows[pos]];
		else F[pos] = Faces[pos_face].Part_Flow[pos_flows[pos]];
		if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) S = Mesh.faces[pos_face].S * Mesh.faces[pos_face].r;
		else S = Mesh.faces[pos_face].S;
		F[pos] = (F[pos] * Mesh.faces[pos_face].n.x + F[pos] * Mesh.faces[pos_face].n.y) * S;
	}
	return F;
}
inline double gm::FISICS::dG(int pos_flow, int iter, bool fluid)
{
	int zone_inlet = 0, zone_outlet = 0, zone_interface = Mesh.bound.interf;
	Mesh.POS_INorOUT(zone_inlet, zone_outlet);
	double G_in = dGi(zone_inlet, pos_flow, fluid), G_out = dGi(zone_outlet, pos_flow, fluid), G_point = 0;
	if (zone_interface != -1) G_point = dGi(zone_interface, pos_flow, fluid);
	if (iter % inf.diter == 0)
	{
		std::cout << INF_FLOW(pos_flow, fluid) << ":   inlet: " << abs(G_in) << ", outlet: " << abs(G_out);
		if (zone_interface != -1) std::cout << ", check_point: " << abs(G_point);
		std::cout << std::endl;
	}
	return abs(abs(G_in) - abs(G_out));
}
std::vector<double> gm::FISICS::dG(int iter, bool fluid)
{
	std::vector<int> pos_flows;
	if (fluid) pos_flows = inf.pos_flows_fluid;
	else pos_flows = inf.pos_flows_part;
	int zone_inlet = 0, zone_outlet = 0, size = (int)pos_flows.size(), zone_interface = Mesh.bound.interf;
	Mesh.POS_INorOUT(zone_inlet, zone_outlet);
	std::vector<double> G_in(size), G_out(size), G_point(size);
	G_in = dGi(zone_inlet, fluid), G_out = dGi(zone_outlet, fluid);
	if (zone_interface != -1) G_point = dGi(zone_interface, fluid);
	if (iter % inf.diter == 0)
	{
		for (int i = 0; i < size; i++)
		{
			std::cout << INF_FLOW(pos_flows[i], fluid) << ":   inlet = " << abs(G_in[i]) << ", outlet = " << abs(G_out[i]);
			if (zone_interface != -1) std::cout << ", check_point = " << abs(G_point[i]);
			std::cout << std::endl;
		}
	}
	for (int i = 0; i < size; i++)
		G_in[i] = abs(abs(G_in[i]) - abs(G_out[i]));
	return G_in;
}
inline double gm::FISICS::Tet()
{
	int zone_inlet = 0, zone_outlet = 0;
	Mesh.POS_INorOUT(zone_inlet, zone_outlet);
	double G_in = dGi(zone_inlet, inf.eq_p, 0), G_out = dGi(zone_outlet, inf.eq_p, 0);
	return 1 - abs(G_out / G_in);
}

void gm::FISICS::FLOW_FACE()
{
	for (int pos_zone = 0; pos_zone < Mesh.zone_faces.size(); pos_zone++) FLOW_FACE(pos_zone);
}

void gm::FISICS::FLOW_FACE(int pos_zone)
{
	switch (Mesh.zone_faces[pos_zone].type)
	{
	case 5: { FLOW_PRESSURE_OUTLET(pos_zone); break; }
	case 14: { FLOW_MASS_FLOW_INLET(pos_zone); break; }
	case 10: { FLOW_VELOCITY_INLET(pos_zone); break; }
	case 24: { FLOW_OUTFLOW(pos_zone); break; }
	default: { FLOW_VAN_LEER(pos_zone); break; }
	}
}

void gm::FISICS::FLOW_PRESSURE_OUTLET(int pos_zone)
{
	int i = Mesh.zone_faces[pos_zone].face_first;
	std::thread F0([&] {
		for (int k = i; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) PRESSURE_OUTLET(k); });
	std::thread F1([&] {
		for (int k = i + 1; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) PRESSURE_OUTLET(k); });
	std::thread F2([&] {
		for (int k = i + 2; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) PRESSURE_OUTLET(k); });
	std::thread F3([&] {
		for (int k = i + 3; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) PRESSURE_OUTLET(k); });
	F0.join(); F1.join(); F2.join(); F3.join();
}

void gm::FISICS::FLOW_MASS_FLOW_INLET(int pos_zone)
{
	int i = Mesh.zone_faces[pos_zone].face_first;
	std::thread F0([&] {
		for (int k = i; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) MASS_FLOW_INLET(k); });
	std::thread F1([&] {
		for (int k = i + 1; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) MASS_FLOW_INLET(k); });
	std::thread F2([&] {
		for (int k = i + 2; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) MASS_FLOW_INLET(k); });
	std::thread F3([&] {
		for (int k = i + 3; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) MASS_FLOW_INLET(k); });
	F0.join(); F1.join(); F2.join(); F3.join();
}

void gm::FISICS::FLOW_VELOCITY_INLET(int pos_zone)
{
	int i = Mesh.zone_faces[pos_zone].face_first;
	std::thread F0([&] {
		for (int k = i; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) VELOCITY_INLET(k); });
	std::thread F1([&] {
		for (int k = i + 1; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) VELOCITY_INLET(k); });
	std::thread F2([&] {
		for (int k = i + 2; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) VELOCITY_INLET(k); });
	std::thread F3([&] {
		for (int k = i + 3; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) VELOCITY_INLET(k); });
	F0.join(); F1.join(); F2.join(); F3.join();
}

void gm::FISICS::FLOW_OUTFLOW(int pos_zone)
{
	int i = Mesh.zone_faces[pos_zone].face_first;
	std::thread F0([&] {
		for (int k = i; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) OUTFLOW(k); });
	std::thread F1([&] {
		for (int k = i + 1; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) OUTFLOW(k); });
	std::thread F2([&] {
		for (int k = i + 2; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) OUTFLOW(k); });
	std::thread F3([&] {
		for (int k = i + 3; k <= Mesh.zone_faces[pos_zone].face_last; k += 4) OUTFLOW(k); });
	F0.join(); F1.join(); F2.join(); F3.join();
}

void gm::FISICS::FLOW_VAN_LEER(int pos_zone)
{
	int i = Mesh.zone_faces[pos_zone].face_first;
	std::thread F0([&] {
		for (int k = i; k <= Mesh.zone_faces[pos_zone].face_last; k += 4)
		{
			VAN_LEER(k);
			if (inf.part) FLOW_PART(k);
		}});
		std::thread F1([&] {
			for (int k = i + 1; k <= Mesh.zone_faces[pos_zone].face_last; k += 4)
			{
				VAN_LEER(k);
				if (inf.part) FLOW_PART(k);
			}});
			std::thread F2([&] {
				for (int k = i + 2; k <= Mesh.zone_faces[pos_zone].face_last; k += 4)
				{
					VAN_LEER(k);
					if (inf.part) FLOW_PART(k);
				}});
				std::thread F3([&] {
					for (int k = i + 3; k <= Mesh.zone_faces[pos_zone].face_last; k += 4)
					{
						VAN_LEER(k);
						if (inf.part) FLOW_PART(k);
					}});
					F0.join(); F1.join(); F2.join(); F3.join();
}

inline void gm::FISICS::PRESSURE_OUTLET(int pos_face)
{
	double U, Ro, P, c, un, k;
	int pos_cell = Mesh.faces[pos_face].cell_l;
	c = val::C(inf, Cells[pos_cell]);
	k = val::k(inf, Cells[pos_cell]);
	un = Un(pos_cell, pos_face);
	U = un, Ro = Cells[pos_cell].Fluid.Opt.Ro, P = Cells[pos_cell].Fluid.Opt.P;
	if (abs(un) < c)
	{
		if (P >= inf.out.P)
		{
			double _C, a;
			a = (k - 1.) / (2. * k) * Cells[pos_cell].Fluid.Opt.Ro * c * (1. - inf.out.P / Cells[pos_cell].Fluid.Opt.P) / (1. - pow(inf.out.P / Cells[pos_cell].Fluid.Opt.P, (k - 1.) / (2. * k)));
			U = un - (inf.out.P - P) / a;
			_C = c + (k - 1.) / 2. * (un - U);
			if (U < _C)
			{
				Ro = k * inf.out.P / (_C * _C);
				P = inf.out.P;
			}
			else
			{
				U = 2. * c / (k + 1.) + (k - 1.) / (k + 1.) * un;
				Ro = Ro * pow(U * U / (c * c), 1. / (k - 1.));
				P = Ro * U * U / k;
			}
		}
		else
		{
			double a;
			a = pow(Ro * (((k + 1.) / 2. * inf.out.P + (k - 1.) / 2. * P)), 0.5);
			U = un - (inf.out.P - P) / a;
			P = inf.out.P;
			if (U >= 0) Ro = Ro * a / (a - Ro * (un - U));
			else Ro = k / (k - 1.) * Cells[pos_cell].Fluid.Opt.P / (k / (k - 1.) * inf.out.P / inf.out.Ro - (un * un + Ut(pos_cell, pos_face) * Ut(pos_cell, pos_face)) / 2.);
		}
	}
	F_FLUID(pos_face, U, Ut(pos_cell, pos_face), P, Ro, 0);
	if (inf.part) Faces[pos_face].Part_Flow = F_PART(pos_cell, pos_face);
}

inline double gm::FISICS::MFI(int pos_face, double Pi)
{
	int pos_cell = Mesh.faces[pos_face].cell_l;
	double u = Un(pos_cell, pos_face), p = Cells[pos_cell].Fluid.Opt.P, U, k = val::k(inf), G;
	U = u - (Pi - p) / pow(Cells[pos_cell].Fluid.Opt.Ro * ((k + 1.) / 2. * Pi + (k - 1.) / 2. * p), 0.5);
	if (inf.Metod_mass_flow == 1) G = inf.in.G;
	else G = inf.in.Rot * inf.in.U01 / inf.in.P1 * pow(Pi, inf.in.nu) * (1 - inf.in.z);
	return k / (k - 1) * Pi * U / G + U * U / 2. - val::Cp(inf) * inf.in.T;
}

void gm::FISICS::MASS_FLOW_INLET(int pos_face)
{
	Meas meas;
	double Ro, U, P, Pi, temp = 0, d = 1, G;
	int pos_cell = Mesh.faces[pos_face].cell_l;
	/*
	c = val::C(inf, Cells[pos_cell]);
	temp = (Cells[pos_cell].Fluid.Opt.Ro * c * Un(pos_cell,pos_face) + Cells[pos_cell].Fluid.Opt.P) / inf.T_0;
	Ro = 0.5 * (temp + pow(temp * temp - 4. * Cells[pos_cell].Fluid.Opt.Ro * c / inf.T_0 * inf.G_0, 0.5));
	U = inf.G_0 / Ro;
	P = Ro * inf.T_0;
	*/
	P = Cells[pos_cell].Fluid.Opt.P;
	for (; d > meas.eps; P = Pi)
	{
		temp = MFI(pos_face, P);
		if (MFI(pos_face, P + meas.eps) - temp == 0) std::cerr << " Error Mass_Flow /= 0", WRITE(), exit(1);
		else Pi = P - (temp * meas.eps) / (MFI(pos_face, P + meas.eps) - temp);
		d = abs(Pi - P);
	}
	Pi = Cells[pos_cell].Fluid.Opt.P;
	U = Un(pos_cell, pos_face) - (P - Pi) / pow(Cells[pos_cell].Fluid.Opt.Ro * (P * (val::k(inf) + 1.) / 2. + Pi * (val::k(inf) - 1.) / 2.), 0.5);
	if (inf.Metod_mass_flow == 1) G = inf.in.G;
	else G = inf.in.Rot * inf.in.U01 / inf.in.P1 * pow(P, inf.in.nu) * (1 - inf.in.z);
	Ro = G / U;

	F_FLUID(pos_face, U, 0., P, Ro);
	if (inf.part) F_PART(pos_face, U, 0., P / (Ro * val::R(inf)), Ro);
}

inline void gm::FISICS::VELOCITY_INLET(int pos_face)
{
	double P = inf.in.P, Ro = inf.in.Ro, U;
	int pos_cell = Mesh.faces[pos_face].cell_l;
	U = inf.in.M * pow(val::k(inf) * P / Ro, 0.5) * Mesh.faces[pos_face].n.x;

	F_FLUID(pos_face, U, 0., P, Ro);
	if (inf.part) F_PART(pos_face, U, 0., P / (Ro * val::R(inf)), Ro);
}

void gm::FISICS::OUTFLOW(int pos_face)
{
	int pos_cell = Mesh.faces[pos_face].cell_l;

	Faces[pos_face].Fluid_Flow = F_FLUID(pos_cell, pos_face);
	T_1(pos_face);
	if (inf.part) Faces[pos_face].Part_Flow = F_PART(pos_cell, pos_face);
}

inline bool gm::FISICS::WALL(int pos_face)
{
	int type = Mesh.TYPE_FACE(pos_face);
	if (type == 3 || type == 7 || type == 25) return 1;
	else return 0;
}

void gm::FISICS::VAN_LEER(int pos_face)
{
	std::vector<double> f_pl(inf.eq_f + inf.eq_fi), f_m = f_pl;
	int pos_cell = Mesh.faces[pos_face].cell_l;
	f_pl = VAN_LEER(pos_cell, pos_face, 1);
	if (!WALL(pos_face)) pos_cell = Mesh.faces[pos_face].cell_r;
	f_m = VAN_LEER(pos_cell, pos_face, -1);
	val::plus(f_pl, f_m);
	Faces[pos_face].Fluid_Flow = f_pl;
	T_1(pos_face);
}

inline std::vector<double> gm::FISICS::VAN_LEER(int pos_cell, int pos_face, int dir)
{
	std::vector<double> f(inf.eq_f + inf.eq_fi);
	double M = val::M(inf, Cells[pos_cell], Un(pos_cell, pos_face));
	if (WALL(pos_face)) M *= dir;
	if (abs(M) < 1)
	{
		double c = val::C(inf, Cells[pos_cell]), k = val::k(inf, Cells[pos_cell]), temp;
		temp = dir * Cells[pos_cell].Fluid.Opt.Ro * c * pow(M + dir, 2) / 4.;
		f[0] = temp;
		f[1] = temp * 2. * c / k * ((k - 1.) / 2. * M + dir);
		f[2] = temp * Ut(pos_cell, pos_face);
		f[3] = temp * (2. * c * c / (k * k - 1.) * pow((k - 1.) / 2. * M + dir, 2) + Ut(pos_cell, pos_face) * Ut(pos_cell, pos_face) / 2.);
		for (int j = 0; j < inf.eq_fi; j++) f[inf.eq_f + j] = f[0] * Cells[pos_cell].Fluid.Opt.Y[j];
		return f;
	}
	if (M * dir >= dir)
		f = F_FLUID(pos_cell, pos_face);
	return f;
}

inline void gm::FISICS::FLOW_PART(int pos_face)
{
	std::vector<double> f(inf.eq_p + inf.eq_pi);
	if (!WALL(pos_face))
	{
		int cell_l = Mesh.faces[pos_face].cell_l, cell_r = Mesh.faces[pos_face].cell_r;
		double Ul = Un(cell_l, pos_face, 0), Ur = Un(cell_r, pos_face, 0);
		if (Ul > 0 && Ur >= 0) f = F_PART(cell_l, pos_face);
		if (Ul <= 0 && Ur < 0) f = F_PART(cell_r, pos_face);
		if (Ul > 0 && Ur < 0)
		{
			double Ro_l = Cells[cell_l].Part.Opt.Ro, Ro_r = Cells[cell_r].Part.Opt.Ro, D;
			if ((pow(Ro_l, 0.5) + pow(Ro_r, 0.5)) == 0) D = 0;
			else D = (pow(Ro_l, 0.5) * Ul + pow(Ro_r, 0.5) * Ur) / (pow(Ro_l, 0.5) + pow(Ro_r, 0.5));
			if (D > 0) f = F_PART(cell_l, pos_face);
			if (D < 0) f = F_PART(cell_r, pos_face);
		}
	}
	Faces[pos_face].Part_Flow = f;
}

void gm::FISICS::FLOW_CELL(double dt)
{
	std::thread F0([&] {
		for (int pos_cell = 0; pos_cell < Cells.size(); pos_cell += 4) FLOW_CELL_MAIN(dt, pos_cell); });
	std::thread F1([&] {
		for (int pos_cell = 1; pos_cell < Cells.size(); pos_cell += 4) FLOW_CELL_MAIN(dt, pos_cell); });
	std::thread F2([&] {
		for (int pos_cell = 2; pos_cell < Cells.size(); pos_cell += 4) FLOW_CELL_MAIN(dt, pos_cell); });
	std::thread F3([&] {
		for (int pos_cell = 3; pos_cell < Cells.size(); pos_cell += 4) FLOW_CELL_MAIN(dt, pos_cell); });
	F0.join(); F1.join(); F2.join(); F3.join();
}

void gm::FISICS::FLOW_CELL_MAIN(double dt, int pos_cell)
{
	FLOW_CELL(dt, pos_cell);
	if (inf.part)
	{
		FLOW_CELL(dt, pos_cell, 0);
		INTERACTION(dt, pos_cell);
		NEW_OPT(pos_cell, 0);
	}
	NEW_OPT(pos_cell);
}

inline void gm::FISICS::FLOW_CELL(double dt, int pos_cell, bool fluid)
{
	double Omega;
	int eq = inf.eq_p + inf.eq_pi;
	if (fluid) eq = inf.eq_f + inf.eq_fi;
	std::vector<double> f(eq), fr;
	if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) Omega = Mesh.cells[pos_cell].Omega_axis;
	else Omega = Mesh.cells[pos_cell].Omega;
	fr = FLOW_CELL(pos_cell, fluid);
	val::mult(fr, -dt / Omega), val::plus(f, fr);
	if (Mesh.bound.axis != -1)
	{
		double r = dt;
		if (inf.Metod_axis == 1) r /= Mesh.cells[pos_cell].r;
		else r *= Mesh.cells[pos_cell].Omega / Omega;
		fr = AXIS(pos_cell, fluid);
		val::mult(fr, r), val::plus(f, fr);
	}
	if (inf.burn)
	{
		double r = dt;
		if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) r *= Mesh.cells[pos_cell].Omega / Omega;
		fr = BURN(pos_cell, fluid);
		val::mult(fr, r), val::plus(f, fr);
	}
	if (fluid) val::plus(Cells[pos_cell].Fluid.Flow, f);
	else val::plus(Cells[pos_cell].Part.Flow, f);
}

std::vector<double> gm::FISICS::FLOW_CELL(int pos_cell, bool fluid)
{
	int eq, pos_face;
	if (fluid) eq = inf.eq_f + inf.eq_fi;
	else eq = inf.eq_p + inf.eq_pi;
	std::vector<double> f(eq), f_temp(eq);
	double S = 0;
	for (int i = 0; i < Mesh.cells[pos_cell].pos_faces.size(); i++)
	{
		pos_face = Mesh.cells[pos_cell].pos_faces[i];
		if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) S = Mesh.faces[pos_face].S * Mesh.faces[pos_face].r;
		else S = Mesh.faces[pos_face].S;
		if (pos_cell == Mesh.faces[pos_face].cell_r) S *= -1;
		if (fluid) f_temp = Faces[pos_face].Fluid_Flow;
		else f_temp = Faces[pos_face].Part_Flow;
		val::mult(f_temp, S), val::plus(f, f_temp);
	}
	return f;
}

inline std::vector<double> gm::FISICS::AXIS(int pos_cell, bool fluid)
{
	std::vector<double> f;
	if (fluid) f = AXIS_FLUID(pos_cell);
	else f = AXIS_PART(pos_cell);
	return f;
}

std::vector<double> gm::FISICS::AXIS_FLUID(int pos_cell)
{
	std::vector<double> f(inf.eq_f + inf.eq_fi);
	switch (inf.Metod_axis)
	{
	case 1:
	{
		f[0] = -Cells[pos_cell].Fluid.Opt.Ro * Cells[pos_cell].Fluid.Opt.V;
		f[1] = -Cells[pos_cell].Fluid.Opt.Ro * Cells[pos_cell].Fluid.Opt.U * Cells[pos_cell].Fluid.Opt.V;
		f[2] = -Cells[pos_cell].Fluid.Opt.Ro * Cells[pos_cell].Fluid.Opt.V * Cells[pos_cell].Fluid.Opt.V;
		f[3] = -Cells[pos_cell].Fluid.Opt.Ro * Cells[pos_cell].Fluid.Opt.V * val::H(inf, Cells[pos_cell]);
		for (int j = 0; j < inf.eq_fi; j++) f[inf.eq_f + j] = f[0] * Cells[pos_cell].Fluid.Opt.Y[j];
		return f;
	}
	default:
	{
		f[2] = Cells[pos_cell].Fluid.Opt.P; return f;
	}
	}
}

std::vector<double> gm::FISICS::AXIS_PART(int pos_cell)
{
	std::vector<double> f(inf.eq_p + inf.eq_pi);
	switch (inf.Metod_axis)
	{
	case 1:
	{
		f[0] = -Cells[pos_cell].Part.Opt.Ns * Cells[pos_cell].Part.Opt.V;
		f[1] = -Cells[pos_cell].Part.Opt.Ro * Cells[pos_cell].Part.Opt.V;
		f[2] = -Cells[pos_cell].Part.Opt.Ro * Cells[pos_cell].Part.Opt.U * Cells[pos_cell].Part.Opt.V;
		f[3] = -Cells[pos_cell].Part.Opt.Ro * Cells[pos_cell].Part.Opt.V * Cells[pos_cell].Part.Opt.V;
		f[4] = -Cells[pos_cell].Part.Opt.Ro * Cells[pos_cell].Part.Opt.V * Cells[pos_cell].Part.Opt.T * val::Cs(inf, Cells[pos_cell]);
		for (int j = 0; j < inf.eq_pi; j++) f[inf.eq_p + j] = f[1] * Cells[pos_cell].Part.Opt.Y[j];
		return f;
	}
	default: return f;
	}
}

std::vector<double> gm::FISICS::BURN(int pos_cell, bool fluid)
{
	int eq = inf.eq_p, eq_i = inf.eq_pi, flow = 1;
	if (fluid) eq = inf.eq_f, eq_i = inf.eq_fi, flow = 0;
	std::vector<double> f(eq + eq_i);
	double Nu_M = 0, Ns_Mdt = Cells[pos_cell].Part.Opt.Ns * val::Mdt(inf, Cells[pos_cell]);
	f[flow] = Ns_Mdt / val::Nu_M(inf, 0, 0);
	for (int i = 0; i < eq_i; i++) f[eq + i] = Ns_Mdt * val::Nu_M(inf, i, fluid) / val::Nu_M(inf, 0, 0), Nu_M += val::Nu_M(inf, i, fluid);
	Nu_M += val::Nu_M(inf, eq_i, fluid);
	f[flow] *= Nu_M;
	if (!fluid) f[eq - 1] = -inf.in.Q * Ns_Mdt;
	return f;
}

inline void gm::FISICS::INTERACTION(double dt, int pos_cell)
{
	std::vector<double> T = Cells[pos_cell].Fluid.Flow, Ts = Cells[pos_cell].Part.Flow;
	double r = dt;
	if (Mesh.bound.axis != -1 && inf.Metod_axis == 2) r *= Mesh.cells[pos_cell].Omega / Mesh.cells[pos_cell].Omega_axis;
	Cells[pos_cell].Fluid.Opt.Ro = T[0];
	if (T[0] <= 0 || T[0] != T[0]) std::cerr << "Error Interection R <= 0 || nan", WRITE(), exit(1);
	if (Ts[0] < 0) Ts[0] = 0;
	if (Ts[1] < 0) Ts[1] = 0;
	Cells[pos_cell].Part.Opt.Ns = Ts[0];
	Cells[pos_cell].Part.Opt.Ro = Ts[1];
	Yi(pos_cell), Yi(pos_cell, 0);
	double dT = Ts[1] / T[0], V1, V2,
		F_1 = val::fi_1(inf, Cells[pos_cell]) * r / (1. + val::fi_1(inf, Cells[pos_cell]) * r),
		F_2 = val::fi_1(inf, Cells[pos_cell]) * r / (1. + dT * val::fi_1(inf, Cells[pos_cell]) * r);
	if (F_1 != F_1) std::cerr << " F_1(fi_1) = nan", WRITE(), exit(1);
	if (F_2 != F_2) std::cerr << " F_2(fi_1) = nan", WRITE(), exit(1);
	//RoU
	Cells[pos_cell].Fluid.Flow[1] = (T[1] + F_1 * Ts[2]) / (1. + dT * F_1);
	Cells[pos_cell].Part.Flow[2] = (Ts[2] + F_2 * T[1] * dT) / (1. + F_2);
	//RoV
	Cells[pos_cell].Fluid.Flow[2] = (T[2] + F_1 * Ts[3]) / (1. + dT * F_1);
	Cells[pos_cell].Part.Flow[3] = (Ts[3] + F_2 * T[2] * dT) / (1. + F_2);
	//RoE or RosTsCs
	dT *= val::Cs(inf, Cells[pos_cell]) / val::Cv(inf, Cells[pos_cell]);
	F_1 = val::fi_2(inf, Cells[pos_cell]) * r / (1. + val::fi_2(inf, Cells[pos_cell]) * r),
		F_2 = val::fi_2(inf, Cells[pos_cell]) * r / (1. + dT * val::fi_2(inf, Cells[pos_cell]) * r);
	if (F_1 != F_1) std::cerr << " F_1(fi_2) = nan", WRITE(), exit(1);
	if (F_2 != F_2) std::cerr << " F_2(fi_2) = nan", WRITE(), exit(1);
	V2 = (pow(Cells[pos_cell].Fluid.Flow[1], 2) + pow(Cells[pos_cell].Fluid.Flow[2], 2)) / (2. * T[0]);
	if (Ts[1] == 0) V1 = 0;
	else
	{
		V1 = val::fi_1(inf, Cells[pos_cell]) * r * (Cells[pos_cell].Part.Flow[2] * (Cells[pos_cell].Fluid.Flow[1] / T[0] - Cells[pos_cell].Part.Flow[2] / Ts[1])
			+ Cells[pos_cell].Part.Flow[3] * (Cells[pos_cell].Fluid.Flow[2] / T[0] - Cells[pos_cell].Part.Flow[3] / Ts[1]));
	}
	Cells[pos_cell].Fluid.Flow[3] = (T[3] + F_1 * Ts[4] - V1 + F_1 * dT * V2) / (1. + F_1 * dT);
	Cells[pos_cell].Part.Flow[4] = (Ts[4] + T[3] * dT * F_2 - dT * F_2 * (V1 + V2)) / (1. + F_2);
}

void gm::FISICS::Yi(int pos_cell, bool fluid)
{
	int i, eq, eq_i;
	double* Y, * Ro_i, Ro, Y_end = 0;
	if (fluid)
	{
		eq = inf.eq_f, eq_i = inf.eq_fi;
		Ro = Cells[pos_cell].Fluid.Flow[0];
		Y = &Cells[pos_cell].Fluid.Opt.Y[0];
		Ro_i = &Cells[pos_cell].Fluid.Flow[eq - 1];
	}
	else
	{
		eq = inf.eq_p, eq_i = inf.eq_pi;
		Ro = Cells[pos_cell].Part.Flow[1];
		Y = &Cells[pos_cell].Part.Opt.Y[0];
		Ro_i = &Cells[pos_cell].Part.Flow[eq - 1];
	}
	for (i = 0; i < eq_i; i++)
	{
		if (Ro == 0 || *(Ro_i + i + 1) < 0) *(Y + i) = 0;
		else *(Y + i) = *(Ro_i + i + 1) / Ro;
		Y_end += *(Y + i);
	}
	if (Ro != 0) *(Y + eq_i) = 1 - Y_end;
	else *(Y + eq_i) = 0;
	if (*(Y + eq_i) < 0) *(Y + eq_i) = 0;
}

inline void gm::FISICS::NEW_OPT(int pos_cell, bool fluid)
{
	if (fluid)
	{
		Cells[pos_cell].Fluid.Opt.Ro = Cells[pos_cell].Fluid.Flow[0];
		if (Cells[pos_cell].Fluid.Opt.Ro <= 0 || Cells[pos_cell].Fluid.Opt.Ro != Cells[pos_cell].Fluid.Opt.Ro)
			std::cerr << " Ro <= 0 || nan ", WRITE(), exit(1);
	}
	else
	{
		Cells[pos_cell].Part.Opt.Ns = Cells[pos_cell].Part.Flow[0];
		if (Cells[pos_cell].Part.Opt.Ns < 0) Cells[pos_cell].Part.Opt.Ns = 0;
		Cells[pos_cell].Part.Opt.Ro = Cells[pos_cell].Part.Flow[1];
		if (Cells[pos_cell].Part.Opt.Ro < 0) Cells[pos_cell].Part.Opt.Ro = 0;
	}
	Yi(pos_cell, fluid);
	NEW_OPT_UV(pos_cell, fluid);
	if (fluid)
	{
		Cells[pos_cell].Fluid.Opt.P = (Cells[pos_cell].Fluid.Flow[3] - (pow(Cells[pos_cell].Fluid.Flow[1], 2) + pow(Cells[pos_cell].Fluid.Flow[2], 2)) / (2. * Cells[pos_cell].Fluid.Flow[0])) * (val::k(inf, Cells[pos_cell]) - 1.);
		if (Cells[pos_cell].Fluid.Opt.P <= 0 || Cells[pos_cell].Fluid.Opt.P != Cells[pos_cell].Fluid.Opt.P)
			std::cerr << " P <= 0 || nan ", WRITE(), exit(1);
	}
	else
	{
		if (Cells[pos_cell].Part.Flow[1] * val::Cs(inf, Cells[pos_cell]) == 0) Cells[pos_cell].Part.Opt.T = 0;
		else Cells[pos_cell].Part.Opt.T = Cells[pos_cell].Part.Flow[4] / (Cells[pos_cell].Part.Flow[1] * val::Cs(inf, Cells[pos_cell]));
	}
}

inline void gm::FISICS::NEW_OPT_UV(int pos_cell, bool fluid)
{
	double* Fl, * u, * v;
	if (fluid)
	{
		Fl = &Cells[pos_cell].Fluid.Flow[0];
		u = &Cells[pos_cell].Fluid.Opt.U, v = &Cells[pos_cell].Fluid.Opt.V;
	}
	else
	{
		Fl = &Cells[pos_cell].Part.Flow[1];
		u = &Cells[pos_cell].Part.Opt.U, v = &Cells[pos_cell].Part.Opt.V;
	}
	if (!fluid && *Fl == 0) *u = 0, * v = 0;
	else *u = *(Fl + 1) / (*Fl), * v = *(Fl + 2) / (*Fl);
}