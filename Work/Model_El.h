#pragma once
#include <vector>
namespace gm
{
	struct NULL_VEC
	{
		double x = 0., y = 0.;
	};

	struct POINT
	{
		int pos = 0;
		double x = 0., y = 0., z = 0.;
		bool type = 0;
	};

	struct FACE
	{
		NULL_VEC centre, n, t;
		int pos = 0, cell_l = 0, cell_r = 0, el_type = 0;
		double S = 0., r = 0.;
		std::vector<int> pos_nodes;
	};

	struct ZONE_FACES
	{
		int face_first = 0, face_last = 0, type = 0;
	};

	struct CELL
	{
		int pos = 0;
		NULL_VEC centre;
		double Omega = 0., Omega_axis = 0., h = 0., r = 0.;
		std::vector<int> pos_faces, pos_nodes;
	};

	struct ZONE_CELLS
	{
		int cell_first = 0, cell_last = 0, type = 0, el_type = 0;
	};

	struct BOUND_TYPE
	{
		int axis = -1, wall = -1, inter = -1, press_out = -1, symm = -1,
			vel_in = -1, mass_flow_in = -1, interf = -1, outflow = -1;
	};
}