#pragma once
#include "Model_El.h"
#include <string>
#include <fstream>

namespace gm
{
	class MODEL_MESH
	{
	public:
		BOUND_TYPE bound;
		int count_all_points, count_all_cells, count_all_faces, dimension;
		std::vector<POINT> points;
		std::vector<FACE> faces;
		std::vector<CELL> cells;
		std::vector<ZONE_FACES> zone_faces;
		std::vector<ZONE_CELLS> zone_cells;
		MODEL_MESH();
	public:
		int TYPE_FACE(int pos_face);
		std::string TYPE_FACES(int pos_zone);
		int TYPE_POS(int type);

		int POS_FACE_IN_CELL(int pos_cell, int pos_face);
		int POS_ZONE_FACES(int pos_face);
		void POS_INorOUT(int& pos_inlet, int& pos_outlet);
	private:
		int TRANSLATE_CHAR(char ch);
		int TRANSLATE_INT(std::string s);
		void INPUT_POINT(std::ifstream& File);
		void INPUT_FACES(std::ifstream& File);
		void INPUT_TYPE_FACES(int type_face, int pos);
		int INPUT_FACE(std::ifstream& File, ZONE_FACES& zone);
		void INPUT_CELLS(std::ifstream& File);

		void STRUCT_FACE();
		void CELL_POINT(int pos_cell);
		void DIRECT_POINTS(int pos_cell);
		void STRUCT();
		void STRUCT_CELL();

		bool NEIGHBOR(int pos_face, int pos_face_n);

		void CENTRE(NULL_VEC& centre, std::vector <int> pos_nodes);
		void CHORD(int pos_cell);
		void PATH(int pos_face);
		void RADIUS(NULL_VEC centre, double& r);
		void SQUARE(int pos_cell);

		NULL_VEC VECTOR(int pos_face);
		NULL_VEC VECTOR(int node_0, NULL_VEC node_1);
		double DET(NULL_VEC vec_0, NULL_VEC vec_1);
		double N(double vec_x, double vec_y, int sign = 1);
		void N_TG(int pos_face);
	};
}

