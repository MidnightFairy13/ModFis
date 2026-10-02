#pragma once
#include <string>
#include "Fisics.h"
struct File
{
	Meas meas;
	std::string address = "E:\\File_Work", mod = "Engine\\1_6.msh", temp = "PCoM",
		add, mesh, fluid, part, gas, met, model;
	int inter = 15;
	bool check = 0;
	File();
	void Interval(std::ofstream& file, int cycle, int mod = 2);
	void Change(std::string name);
	void Mod_Circuit(gm::MODEL_MESH& Model);
	void Inf_Cell(gm::MODEL_MESH& Model);
	void Model(gm::MODEL_MESH& Model);
	void Opt_Fluid(std::string name, gm::FISICS& Fis, int pos_cell, int cycle);
	void Opt_Fluid_Flow(std::string name, gm::FISICS& Fis, int pos_cell);
	void Opt_Part(std::string name, gm::FISICS& Fis, int pos_cell, int cycle);
	void Opt_Part_Flow(std::string name, gm::FISICS& Fis, int pos_cell);
	void Opt_Gas(std::string name, gm::FISICS& Fis, int pos_cell, int pos_mat, int cycle);
	void Opt_Metal(std::string name, gm::FISICS& Fis, int pos_cell, int pos_mat, int cycle);
	void Opt_Cell(gm::FISICS& Fis);
	void Opt_Face(gm::FISICS& Fis, int pos_zone);
	void Opt_Flow(gm::FISICS& Fis);
	void Inf_Fis(gm::FISICS& Fis);
	void Inf_Conditions(gm::FISICS& Fis);
	void Opt(gm::FISICS& Fis);
	bool Clear_File(std::string file);
	void Create_File();
};