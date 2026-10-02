#pragma once
#include "File.h"
#include <iostream>
#include <filesystem>
File::File()
{
	add = address + "\\" + temp, mesh = address + "\\Database\\Mesh\\" + mod,
		fluid = add + "\\Fluid", part = add + "\\Part",
		gas = add + "\\Gas", met = add + "\\Metal";
}
void File::Interval(std::ofstream& file, int cycle, int mod)
{
	int size = file.tellp(), dl = inter - (size - cycle * mod) % inter;
	for (int i = 0; i < dl; i++) file << ' ';
}
inline void File::Change(std::string name)
{
	temp = name,
		add = address + "\\" + temp,
		fluid = add + "\\Fluid", part = add + "\\Part",
		gas = add + "\\Gas", met = add + "\\Metal";
}
void File::Mod_Circuit(gm::MODEL_MESH& Model)
{
	std::ofstream file(add + "\\Circuit_Model.txt");
	int point = 0;
	for (int pos_point = 0; pos_point < Model.points.size(); pos_point++)
	{
		if (!Model.points[pos_point].type)
		{
			file << Model.points[pos_point].x * meas.X, Interval(file, point);
			file << Model.points[pos_point].y * meas.X, Interval(file, point);
			file << std::endl;
			point++;
		}
	}
	file << point;
	file.close();
}
void File::Inf_Cell(gm::MODEL_MESH& Model)
{
	std::ofstream inf(add + "\\Inf_Cells.txt");
	for (int pos_cell = 0; pos_cell < Model.cells.size(); pos_cell++)
	{
		inf << Model.cells[pos_cell].pos + 1 << " cell (" << Model.cells[pos_cell].centre.x << ", " << Model.cells[pos_cell].centre.y << ")" << std::endl
			<< "S = " << Model.cells[pos_cell].Omega << ", r = " << Model.cells[pos_cell].r << std::endl;
		if (Model.bound.axis != -1) inf << "S(axis) = " << Model.cells[pos_cell].Omega_axis << ", r = " << Model.cells[pos_cell].r << std::endl;
		inf << "faces: " << std::endl;
		for (int i = 0; i < Model.cells[pos_cell].pos_faces.size(); i++)
		{
			int pos_face = Model.cells[pos_cell].pos_faces[i];
			inf << pos_face + 1 << "- points ( ";
			for (int j = 0; j < Model.faces[pos_face].pos_nodes.size(); j++)
				inf << Model.faces[pos_face].pos_nodes[j] + 1 << " ";
			inf << "), n(" << Model.faces[pos_face].n.x << ", " << Model.faces[pos_face].n.y << "), t(" << Model.faces[pos_face].t.x << ", " << Model.faces[pos_face].t.y << ")" << ", l = " << Model.faces[pos_face].S << ", r = " << Model.faces[pos_face].r << "; cell_l = " << Model.faces[pos_face].cell_l + 1 << ", cell_r = " << Model.faces[pos_face].cell_r + 1 << std::endl;
		}
		inf << "nodes: ";
		for (int i = 0; i < Model.cells[pos_cell].pos_nodes.size(); i++) inf << Model.cells[pos_cell].pos_nodes[i] + 1 << " ";
		inf << std::endl;
	}
	inf.close();
}
inline void File::Model(gm::MODEL_MESH& Model)
{
	Mod_Circuit(Model);
	Inf_Cell(Model);
}
inline void File::Opt_Fluid(std::string name, gm::FISICS& Fis, int pos_cell, int cycle)
{
	std::ofstream file(fluid + "\\" + name + "_(x_y_U_V_P_Ro_T_M).txt", std::ios::app);
	file << Fis.Mesh.cells[pos_cell].centre.x * meas.X, Interval(file, cycle);
	file << Fis.Mesh.cells[pos_cell].centre.y * meas.X, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.U * meas.V, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.V * meas.V, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.P * meas.P, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.Ro * meas.Ro, Interval(file, cycle);
	file << val::T(Fis.inf, Fis.Cells[pos_cell]) * meas.T, Interval(file, cycle);
	file << val::M(Fis.inf, Fis.Cells[pos_cell]), Interval(file, cycle);
	file << std::endl;
	file.close();
}
void File::Opt_Fluid_Flow(std::string name, gm::FISICS& Fis, int pos_cell)
{
	std::ofstream file(fluid + "\\" + name + "_Flow.txt", std::ios::app);
	double* cell, * face;
	cell = &Fis.Cells[pos_cell].Fluid.Flow[0];
	file << std::endl << "Cell " << pos_cell + 1 << ':' << std::endl;
	for (int i = 0; i < Fis.inf.eq_f + Fis.inf.eq_fi; i++) file << *(cell + i) << std::endl;
	for (int i = 0; i < Fis.Mesh.cells[pos_cell].pos_faces.size(); i++)
	{
		int pos_face = Fis.Mesh.cells[pos_cell].pos_faces[i];
		file << std::endl << "     " << "Face " << pos_face + 1 << ':' << std::endl;
		face = &Fis.Faces[pos_face].Fluid_Flow[0];
		for (int j = 0; j < Fis.inf.eq_f + Fis.inf.eq_fi; j++) file << "     " << *(face + j) << std::endl;
	}
}
inline void File::Opt_Part(std::string name, gm::FISICS& Fis, int pos_cell, int cycle)
{
	std::ofstream file(part + "\\" + name + "_(x_y_U_V_Ro_T_Ns_D).txt", std::ios::app);
	file << Fis.Mesh.cells[pos_cell].centre.x * meas.X, Interval(file, cycle);
	file << Fis.Mesh.cells[pos_cell].centre.y * meas.X, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.U * meas.V, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.V * meas.V, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.Ro * meas.Ro, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.T * meas.T, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.Ns / (meas.X * meas.X * meas.X), Interval(file, cycle);
	file << val::Ds(Fis.inf, Fis.Cells[pos_cell]) * meas.D, Interval(file, cycle);
	file << std::endl;
	file.close();
}
void File::Opt_Part_Flow(std::string name, gm::FISICS& Fis, int pos_cell)
{
	std::ofstream file(part + "\\" + name + "_Flow.txt", std::ios::app);
	double* cell, * face;
	cell = &Fis.Cells[pos_cell].Part.Flow[0];
	file << std::endl << "Cell " << pos_cell + 1 << ':' << std::endl;
	for (int i = 0; i < Fis.inf.eq_p + Fis.inf.eq_pi; i++) file << *(cell + i) << std::endl;
	for (int i = 0; i < Fis.Mesh.cells[pos_cell].pos_faces.size(); i++)
	{
		int pos_face = Fis.Mesh.cells[pos_cell].pos_faces[i];
		file << std::endl << "     " << "Face " << pos_face + 1 << ':' << std::endl;
		face = &Fis.Faces[pos_face].Part_Flow[0];
		for (int j = 0; j < Fis.inf.eq_p + Fis.inf.eq_pi; j++) file << "     " << *(face + j) << std::endl;
	}
}
inline void File::Opt_Gas(std::string name, gm::FISICS& Fis, int pos_cell, int pos_mat, int cycle)
{
	std::ofstream file(gas + "\\" + name + '_' + Fis.inf.GS[pos_mat].name + "_(x_y_Ro_Y_X).txt", std::ios::app);
	file << Fis.Mesh.cells[pos_cell].centre.x * meas.X, Interval(file, cycle);
	file << Fis.Mesh.cells[pos_cell].centre.y * meas.X, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.Y[pos_mat] * Fis.Cells[pos_cell].Fluid.Opt.Ro, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Fluid.Opt.Y[pos_mat], Interval(file, cycle);
	file << val::C_mat(Fis.inf, Fis.Cells[pos_cell], pos_mat), Interval(file, cycle);
	file << std::endl;
	file.close();
}

inline void File::Opt_Metal(std::string name, gm::FISICS& Fis, int pos_cell, int pos_mat, int cycle)
{
	std::ofstream file(met + "\\" + name + '_' + Fis.inf.MT[pos_mat].name + "_(x_y_Ro_Y_D).txt", std::ios::app);
	file << Fis.Mesh.cells[pos_cell].centre.x * meas.X, Interval(file, cycle);
	file << Fis.Mesh.cells[pos_cell].centre.y * meas.X, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.Y[pos_mat] * Fis.Cells[pos_cell].Part.Opt.Ro, Interval(file, cycle);
	file << Fis.Cells[pos_cell].Part.Opt.Y[pos_mat], Interval(file, cycle);
	file << val::Dsi(Fis.inf, Fis.Cells[pos_cell], pos_mat) * meas.D, Interval(file, cycle);
	file << std::endl;
	file.close();
}

void File::Opt_Cell(gm::FISICS& Fis)
{
	std::string name = "Cell";
	for (int pos_cell = 0; pos_cell < Fis.Cells.size(); pos_cell++)
	{
		Opt_Fluid(name, Fis, pos_cell, pos_cell);
		for (int pos_gas = 0; pos_gas < Fis.inf.GS.size(); pos_gas++) Opt_Gas(name, Fis, pos_cell, pos_gas, pos_cell);
		if (Fis.inf.part)
		{
			Opt_Part(name, Fis, pos_cell, pos_cell);
			for (int pos_met = 0; pos_met < Fis.inf.MT.size(); pos_met++) Opt_Metal(name, Fis, pos_cell, pos_met, pos_cell);
		}
	}
}
void File::Opt_Face(gm::FISICS& Fis, int pos_zone)
{
	std::string name = Fis.Mesh.TYPE_FACES(pos_zone);
	for (int pos_face = Fis.Mesh.zone_faces[pos_zone].face_first, cycle = 0; pos_face <= Fis.Mesh.zone_faces[pos_zone].face_last; pos_face++, cycle++)
	{
		int pos_cell = Fis.Mesh.faces[pos_face].cell_l;
		Opt_Fluid(name, Fis, pos_cell, cycle);
		for (int pos_gas = 0; pos_gas < Fis.inf.GS.size(); pos_gas++) Opt_Gas(name, Fis, pos_cell, pos_gas, cycle);
		if (Fis.inf.part)
		{
			Opt_Part(name, Fis, pos_cell, cycle);
			for (int pos_met = 0; pos_met < Fis.inf.MT.size(); pos_met++) Opt_Metal(name, Fis, pos_cell, pos_met, cycle);
		}
	}
}

void File::Opt_Flow(gm::FISICS& Fis)
{
	std::string name = "Cell";
	for (int pos_cell = 0; pos_cell < Fis.Cells.size(); pos_cell++)
	{
		Opt_Fluid_Flow(name, Fis, pos_cell);
		if (Fis.inf.part) Opt_Part_Flow(name, Fis, pos_cell);
	}
}

void File::Inf_Fis(gm::FISICS& Fis)
{
	std::ofstream file(add + "\\" + "Inf_flow.txt");
	file << "Flow fluid" << std::endl;
	int inlet, outlet;
	Fis.Mesh.POS_INorOUT(inlet, outlet);
	for (int i = 0; i < Fis.inf.eq_f + Fis.inf.eq_fi; i++)
		file << Fis.INF_FLOW(i) << ": inlet = " << abs(Fis.dGi(inlet, i)) << ", outlet = " << abs(Fis.dGi(outlet, i)) << std::endl;
	if (Fis.inf.part)
	{
		file << std::endl << "Flow part" << std::endl;
		for (int i = 0; i < Fis.inf.eq_p + Fis.inf.eq_pi; i++)
			file << Fis.INF_FLOW(i, 0) << ": inlet = " << abs(Fis.dGi(inlet, i, 0)) << ", outlet = " << abs(Fis.dGi(outlet, i, 0)) << std::endl;
		if (Fis.inf.burn) file << std::endl << "Tet = " << Fis.Tet() << std::endl;
	}
	file.close();
}

inline void File::Inf_Conditions(gm::FISICS& Fis)
{
	std::ofstream file(add + "\\" + "Inf_Conditions.txt");
	file << mesh << std::endl << std::endl;
	int inlet, outlet;
	Fis.Mesh.POS_INorOUT(inlet, outlet);
	file << "dt: " << Fis.inf.ddt << std::endl;
	file << "part: " << Fis.inf.part << std::endl;
	file << "burn: " << Fis.inf.burn << std::endl;
	if (Fis.Mesh.bound.axis != -1) file << "Metod_axis: " << Fis.inf.Metod_axis << std::endl;
	if (Fis.Mesh.bound.mass_flow_in != -1) file << "Metod_mass_flow_inlet: " << Fis.inf.Metod_mass_flow << std::endl;
	if (Fis.inf.burn) file << "Metod_Mdt: " << Fis.inf.Metod_Mdt << std::endl;
	file << std::endl << Fis.Mesh.TYPE_FACES(inlet);
	if (Fis.Mesh.bound.mass_flow_in != -1)
	{
		if (Fis.inf.Metod_mass_flow == 1)
		{
			file << " (G, T): " << std::endl;
			file << Fis.inf.in.G * meas.G << "   " << Fis.inf.in.T * meas.T << std::endl;
		}
		else
		{
			file << " (Rot, U01, nu, Tp): " << std::endl;
			file << Fis.inf.in.Rot * meas.Ro << "   " << Fis.inf.in.U01 * meas.V << "   " << Fis.inf.in.nu << "   " << Fis.inf.in.T * meas.T << std::endl;
		}
	}
	if (Fis.Mesh.bound.vel_in != -1)
	{
		file << " (P, M, Ro): " << std::endl;
		file << Fis.inf.in.P * meas.P << "   " << Fis.inf.in.M << "   " << Fis.inf.in.Ro * meas.Ro << std::endl;
	}
	if (Fis.inf.part)
	{
		file << " (D, z";
		if (Fis.inf.burn) file << ", Q";
		file << "): " << std::endl;
		file << Fis.inf.in.D * meas.D << "   " << Fis.inf.in.z;
		if (Fis.inf.burn) file << "   " << Fis.inf.in.Q * meas.Q;
		file << std::endl;
	}
	file << std::endl << Fis.Mesh.TYPE_FACES(outlet);
	if (Fis.Mesh.bound.press_out != -1)
	{
		file << " (P, Ro): " << std::endl;
		file << Fis.inf.out.P * meas.P << "   " << Fis.inf.out.Ro * meas.Ro << std::endl;
	}
	file << std::endl << "Fluid_Components: " << std::endl;
	for (int i = 0; i < Fis.inf.GS.size(); i++)
		file << Fis.inf.GS[i].name << " (reag, Nu, Bound_Y): " << Fis.inf.GS[i].reag << "   " << Fis.inf.GS[i].Nu << "   " << Fis.inf.GS[i].Bound_Y << std::endl;
	if (Fis.inf.part)
	{
		file << std::endl << "Part_Components: " << std::endl;
		for (int i = 0; i < Fis.inf.MT.size(); i++)
			file << Fis.inf.MT[i].name << " (reag, Nu, Bound_Y): " << Fis.inf.MT[i].reag << "   " << Fis.inf.MT[i].Nu << "   " << Fis.inf.MT[i].Bound_Y << std::endl;
	}
	file.close();
}

void File::Opt(gm::FISICS& Fis)
{
	Create_File();
	std::cout << "Writing Model File...";
	Model(Fis.Mesh); std::cout << "OK" << std::endl;
	std::cout << "Writing Inform File...";
	Inf_Fis(Fis);
	Inf_Conditions(Fis);
	std::cout << "OK" << std::endl;
	std::cout << "Writing Options Cell File...";
	Opt_Cell(Fis); std::cout << "OK" << std::endl;
	std::cout << "Writing Options Face File...";
	for (int pos_zone = 0; pos_zone < Fis.Mesh.zone_faces.size(); pos_zone++)
		if (Fis.Mesh.zone_faces[pos_zone].type != 2 && Fis.Mesh.zone_faces[pos_zone].type != 18) Opt_Face(Fis, pos_zone);
	std::cout << "OK" << std::endl;
	if (check) Opt_Flow(Fis);
}

inline bool File::Clear_File(std::string file)
{
	std::filesystem::path File = file;
	if (std::filesystem::exists(File) && std::filesystem::is_directory(File))
	{
		bool del = 0;
		std::cout << file << ": Delete? (0 or 1) ", std::cin >> del;
		if (del) std::filesystem::remove_all(File), std::cout << "Directory delete" << std::endl;
		else return 1;
	}
	return 0;
}

inline void File::Create_File()
{
	std::filesystem::path File = add;
	bool del = 0;
	del = Clear_File(add);
	if (!del)
	{
		std::cout << "Writing..." << std::endl;
		File = fluid, std::filesystem::create_directories(File);
		File = gas, std::filesystem::create_directory(File);
		File = part, std::filesystem::create_directory(File);
		File = met, std::filesystem::create_directory(File);
	}
	else
	{
		std::string name;
		std::cout << "Input name Directory: ", std::cin >> name;
		File = name;
		if (std::filesystem::exists(File) && std::filesystem::is_directory(File)) std::cout << " Directory Exists!";
		else Change(name);
		Create_File();
	}
}
