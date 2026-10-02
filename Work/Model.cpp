#pragma once
#include "Model.h"
#include "Null_Cond.h"
#include "File.h"
#include <iostream>
gm::MODEL_MESH::MODEL_MESH()
{
	File f;
	std::ifstream file(f.mesh);
	if (!file) { std::cerr << " FILE NOT FOUND!"; file.close(); exit(1); }
	std::string s;
	bound.axis = -1, bound.wall = -1, bound.inter = -1, bound.press_out = -1,
		bound.symm = -1, bound.vel_in = -1, bound.mass_flow_in = -1, bound.interf = -1, bound.outflow = -1;
	std::cout << " Reading Model File...";
	getline(file, s);
	//размерность
	for (; file && s != "(0 \"Dimension:\")";) getline(file, s);
	file >> s, file >> s, s = s[0];
	dimension = TRANSLATE_INT(s);
	//точки
	for (; file && s != "(10";) file >> s; getline(file, s);
	INPUT_POINT(file);
	//грани
	for (; file && s != "(0 \"Faces:\")";) getline(file, s);
	INPUT_FACES(file);
	//ячейки
	for (; file && s != "(0 \"Cells:\")";) getline(file, s);
	INPUT_CELLS(file);
	file.close();
	STRUCT();
	std::cout << " OK" << std::endl;
}

inline int gm::MODEL_MESH::TRANSLATE_CHAR(char ch)
{
	int c;
	if (ch == 'a' || ch == 'A') return 10;
	if (ch == 'b' || ch == 'B') return 11;
	if (ch == 'c' || ch == 'C') return 12;
	if (ch == 'd' || ch == 'D') return 13;
	if (ch == 'e' || ch == 'E') return 14;
	if (ch == 'f' || ch == 'F') return 15;
	c = (int)(ch - '0');
	return c;
}

int gm::MODEL_MESH::TRANSLATE_INT(std::string _s)
{
	int count = 0;
	for (int i = 0; _s[i]; i++) count += TRANSLATE_CHAR(_s[_s.size() - i - 1]) * (int)pow(16, i);
	return count;
}

std::string gm::MODEL_MESH::TYPE_FACES(int pos_zone)
{
	int type = zone_faces[pos_zone].type;
	std::string s;
	switch (type)
	{
	case 2: return "Interior";
	case 3: return "Wall";
	case 5: return "Pressure-outlet";
	case 7: return "Symmetry";
	case 10: return "Velocity-inlet";
	case 14: return "Mass-flow-inlet";
	case 18: return "Interface";
	case 24: return "Outflow";
	case 25: return "Axis";
	default: return "Don't found type faces";
	}
}

int gm::MODEL_MESH::TYPE_POS(int type)
{
	switch (type)
	{
	case 2: return bound.inter;
	case 3: return bound.wall;
	case 5: return bound.press_out;
	case 7: return bound.symm;
	case 10: return bound.vel_in;
	case 14: return bound.mass_flow_in;
	case 18: return bound.interf;
	case 24: return bound.outflow;
	case 25: return bound.axis;
	default: { std::cout << "Error type face", exit(1); }
	}
}

int gm::MODEL_MESH::TYPE_FACE(int pos_face)
{
	int i = POS_ZONE_FACES(pos_face);
	return zone_faces[i].type;
}

int gm::MODEL_MESH::POS_FACE_IN_CELL(int pos_cell, int pos_face)
{
	for (int i = 0; i < cells[pos_cell].pos_faces.size(); i++)
	{
		if (cells[pos_cell].pos_faces[i] == pos_face) return i;
	}
	std::cerr << " Error, don't found index cell" << std::endl, exit(1);
}

int gm::MODEL_MESH::POS_ZONE_FACES(int pos_face)
{
	if (pos_face < 0 && pos_face >= count_all_faces) std::cerr << " Error, don't found index face" << std::endl, exit(1);
	for (int i = 0; i < zone_faces.size(); i++)
		if (pos_face >= zone_faces[i].face_first && pos_face <= zone_faces[i].face_last) return i;
}

void gm::MODEL_MESH::POS_INorOUT(int& pos_inlet, int& pos_outlet)
{
	int pos_zone = 0;
	if (bound.mass_flow_in != -1 || bound.vel_in != -1)
	{
		if (bound.mass_flow_in != -1) pos_zone = bound.mass_flow_in;
		else pos_zone = bound.vel_in;
	}
	else std::cerr << "Error inlet faces", exit(1);
	pos_inlet = pos_zone;
	if (bound.outflow != -1 || bound.press_out != -1)
	{
		if (bound.outflow != -1) pos_zone = bound.outflow;
		else pos_zone = bound.press_out;
	}
	else std::cerr << "Error outlet faces", exit(1);
	pos_outlet = pos_zone;
}

bool gm::MODEL_MESH::NEIGHBOR(int face, int face_n)
{
	if (faces[face].pos_nodes[0] == faces[face_n].pos_nodes[0] || faces[face].pos_nodes[0] == faces[face_n].pos_nodes[1]
		|| faces[face].pos_nodes[1] == faces[face_n].pos_nodes[0] || faces[face].pos_nodes[1] == faces[face_n].pos_nodes[1]) return 1;
	else return 0;
}

void gm::MODEL_MESH::CENTRE(NULL_VEC& centre, std::vector <int> pos_nodes)
{
	double x, y;
	int i;
	for (x = 0, y = 0, i = 0; i < pos_nodes.size(); i++)
	{
		x += points[pos_nodes[i]].x;
		y += points[pos_nodes[i]].y;
	}
	centre.x = x / pos_nodes.size();
	centre.y = y / pos_nodes.size();
}

void gm::MODEL_MESH::CHORD(int pos_cell)
{
	int i, j, pos_face, pos_point_k, pos_point_0, pos_point_1;
	double h = 100, h_temp;
	for (i = 0; i < cells[pos_cell].pos_nodes.size(); i++)
	{
		pos_point_k = cells[pos_cell].pos_nodes[i];
		for (j = 0; j < cells[pos_cell].pos_faces.size(); j++)
		{
			pos_face = cells[pos_cell].pos_faces[j];
			pos_point_0 = faces[pos_face].pos_nodes[0];
			pos_point_1 = faces[pos_face].pos_nodes[1];
			if (pos_point_k != pos_point_0 && pos_point_k != pos_point_1)
			{
				double a = (points[pos_point_1].y - points[pos_point_0].y) / (points[pos_point_0].x - points[pos_point_1].x);
				h_temp = abs(points[pos_point_k].y - points[pos_point_0].y + a * (points[pos_point_k].x - points[pos_point_0].x)) / pow(1. + a * a, 1 / 2);
				if (h > h_temp) h = h_temp;
			}
		}
	}
	cells[pos_cell].h = h;
}

void gm::MODEL_MESH::PATH(int pos_face)
{
	double dx, dy;
	dx = points[faces[pos_face].pos_nodes[0]].x - points[faces[pos_face].pos_nodes[1]].x;
	dy = points[faces[pos_face].pos_nodes[0]].y - points[faces[pos_face].pos_nodes[1]].y;
	faces[pos_face].S = pow(dx * dx + dy * dy, 0.5);
}

void gm::MODEL_MESH::RADIUS(NULL_VEC centre, double& r)
{
	if (bound.axis == -1) r = 1;
	else
	{
		int pos_axis = zone_faces[bound.axis].face_first;
		if (points[faces[pos_axis].pos_nodes[0]].y - points[faces[pos_axis].pos_nodes[1]].y != 0) std::cerr << "Don't correct axis", exit(1);
		r = abs(points[faces[pos_axis].pos_nodes[0]].y - centre.y);
	}
}

void gm::MODEL_MESH::SQUARE(int pos_cell)
{
	if (cells[pos_cell].pos_faces.size() == 3) cells[pos_cell].Omega = 0.5 * abs(DET(VECTOR(cells[pos_cell].pos_faces[0]), VECTOR(cells[pos_cell].pos_faces[1])));
	else
	{
		double S1, S2;
		if (NEIGHBOR(cells[pos_cell].pos_faces[0], cells[pos_cell].pos_faces[1]))
		{
			S1 = 0.5 * abs(DET(VECTOR(cells[pos_cell].pos_faces[0]), VECTOR(cells[pos_cell].pos_faces[1])));
			S2 = 0.5 * abs(DET(VECTOR(cells[pos_cell].pos_faces[2]), VECTOR(cells[pos_cell].pos_faces[3])));
		}
		else
		{
			S1 = 0.5 * abs(DET(VECTOR(cells[pos_cell].pos_faces[0]), VECTOR(cells[pos_cell].pos_faces[2])));
			S2 = 0.5 * abs(DET(VECTOR(cells[pos_cell].pos_faces[1]), VECTOR(cells[pos_cell].pos_faces[3])));
		}
		cells[pos_cell].Omega = S1 + S2;
	}
	if (bound.axis != -1)
	{
		Meas meas;
		std::vector<double> r, h;
		for (int pos_node = 0; pos_node < cells[pos_cell].pos_nodes.size(); pos_node++)
			h.push_back(points[cells[pos_cell].pos_nodes[pos_node]].x), r.push_back(points[cells[pos_cell].pos_nodes[pos_node]].y);
		if (cells[pos_cell].pos_nodes.size() == 3) h.push_back(h[0]), r.push_back(r[0]);
		cells[pos_cell].Omega_axis = meas.pi / 3. * abs(((r[3] - r[1]) * (h[2] * (r[1] + r[3] + r[2]) - h[0] * (r[1] + r[3] + r[0])) - (r[2] - r[0]) * (h[3] * (r[0] + r[3] + r[2]) - h[1] * (r[0] + r[1] + r[2]))));
	}
	else cells[pos_cell].Omega_axis = 0;
}

gm::NULL_VEC gm::MODEL_MESH::VECTOR(int pos_face)
{
	if (faces[pos_face].el_type != 2) std::cerr << " Error, el_type != 2" << std::endl, exit(1);
	NULL_VEC vec;
	vec.x = points[faces[pos_face].pos_nodes[1]].x - points[faces[pos_face].pos_nodes[0]].x;
	vec.y = points[faces[pos_face].pos_nodes[1]].y - points[faces[pos_face].pos_nodes[0]].y;
	return vec;
}

gm::NULL_VEC gm::MODEL_MESH::VECTOR(int pos_node_0, NULL_VEC node_1)
{
	NULL_VEC vec;
	vec.x = node_1.x - points[pos_node_0].x;
	vec.y = node_1.y - points[pos_node_0].y;
	return vec;
}

double gm::MODEL_MESH::DET(NULL_VEC vec_0, NULL_VEC vec_1)
{
	return vec_0.x * vec_1.y - vec_0.y * vec_1.x;
}

double gm::MODEL_MESH::N(double vec_x, double vec_y, int sign)
{
	if (vec_y == 0) return 0.;
	return sign * vec_y / pow(vec_y * vec_y + vec_x * vec_x, 0.5);
}

void gm::MODEL_MESH::N_TG(int pos_face)
{
	NULL_VEC edge;
	edge = VECTOR(pos_face);
	faces[pos_face].n.x = N(edge.x, edge.y, 1);
	faces[pos_face].t.y = faces[pos_face].n.x;
	faces[pos_face].n.y = N(edge.y, edge.x, -1);
	faces[pos_face].t.x = -faces[pos_face].n.y;
}

void gm::MODEL_MESH::INPUT_POINT(std::ifstream& File)
{
	Meas meas;
	std::string s;
	int i, j;
	POINT point;
	File >> s, File >> s;
	File >> s, i = TRANSLATE_INT(s) - 1;
	File >> s, j = TRANSLATE_INT(s) - 1;
	count_all_points = j + 1;
	getline(File, s);
	for (; i <= j; i++)
	{
		point.pos = i;
		point.type = 1;
		File >> s, point.x = stod(s) / meas.X;
		File >> s, point.y = stod(s) / meas.X;
		if (dimension == 3) File >> s, point.z = stod(s) / meas.X;
		else point.z = 0;
		points.push_back(point);
	}
}

void gm::MODEL_MESH::INPUT_FACES(std::ifstream& File)
{
	std::string s;
	int i, j, k;
	ZONE_FACES zone;
	File >> s;
	File >> s, i = TRANSLATE_INT(s) - 1;
	File >> s, j = TRANSLATE_INT(s) - 1;
	count_all_faces = j + 1;
	for (k = 0; i <= j; k++)
	{
		File >> s, File >> s;
		File >> s, zone.face_first = TRANSLATE_INT(s) - 1;
		File >> s, zone.face_last = TRANSLATE_INT(s) - 1;
		File >> s, zone.type = stoi(s);
		INPUT_TYPE_FACES(zone.type, k);
		i = INPUT_FACE(File, zone);
		zone_faces.push_back(zone);
	}
}

void gm::MODEL_MESH::INPUT_TYPE_FACES(int type_face, int pos)
{
	switch (type_face)
	{
	case 2: { bound.inter = pos; break; }
	case 3: { bound.wall = pos; break; }
	case 5: { bound.press_out = pos; break; }
	case 7: { bound.symm = pos; break; }
	case 10: { bound.vel_in = pos; break; }
	case 14: { bound.mass_flow_in = pos; break; }
	case 18: { bound.interf = pos; break; }
	case 24: { bound.outflow = pos; break; }
	case 25: { bound.axis = pos; break; }
	default: std::cerr << "Error type_face not found", exit(1);
	}
}

int gm::MODEL_MESH::INPUT_FACE(std::ifstream& File, ZONE_FACES& zone)
{
	std::string s;
	int i, j;
	FACE face;
	File >> s;
	for (i = zone.face_first; i <= zone.face_last; i++)
	{
		face.pos = i;
		File >> s, face.el_type = stoi(s);
		for (j = 0; j < face.el_type; j++)
		{
			int node;
			File >> s, node = TRANSLATE_INT(s) - 1;
			face.pos_nodes.push_back(node);
			if (zone.type != 2 && zone.type != 18) points[node].type = 0;
		}
		File >> s, face.cell_l = TRANSLATE_INT(s) - 1;
		File >> s, face.cell_r = TRANSLATE_INT(s) - 1;
		if (face.cell_l == -1)
		{
			int point;
			face.cell_l = face.cell_r, face.cell_r = -1;
			if (dimension != 2) std::cerr << " Error dimension != 2", exit(1);
			point = face.pos_nodes[0];
			face.pos_nodes[0] = face.pos_nodes[1];
			face.pos_nodes[1] = point;
		}
		faces.push_back(face);
		face.pos_nodes.clear();
	}
	return i;
}

void gm::MODEL_MESH::INPUT_CELLS(std::ifstream& File)
{
	ZONE_CELLS zone;
	zone.cell_last = -1;
	std::string s;
	int i, j;
	File >> s, File >> s;
	File >> s, i = TRANSLATE_INT(s) - 1;
	File >> s, j = TRANSLATE_INT(s) - 1;
	count_all_cells = j + 1;
	File >> s;
	for (; zone.cell_last < j;)
	{
		File >> s, File >> s;
		File >> s, zone.cell_first = TRANSLATE_INT(s) - 1;
		File >> s, zone.cell_last = TRANSLATE_INT(s) - 1;
		File >> s, zone.type = stoi(s);
		File >> s[0], zone.el_type = stoi(s);
		File >> s;
		for (i = zone.cell_first; i <= zone.cell_last; i++)
		{
			CELL cell;
			cell.pos = i;
			cells.push_back(cell);
		}
		zone_cells.push_back(zone);
	}
}

void gm::MODEL_MESH::STRUCT_FACE()
{
	int cl, cr;
	for (int pos_face = 0; pos_face < count_all_faces; pos_face++)
	{
		cl = faces[pos_face].cell_l;
		if (cl != -1) cells[cl].pos_faces.push_back(pos_face);
		else std::cerr << "Error cl = -1", exit(1);
		cr = faces[pos_face].cell_r;
		if (cr != -1) cells[cr].pos_faces.push_back(pos_face);
		if (dimension == 2)
		{
			PATH(pos_face);
			CENTRE(faces[pos_face].centre, faces[pos_face].pos_nodes);
			RADIUS(faces[pos_face].centre, faces[pos_face].r);
			N_TG(pos_face);
		}
	}
}

void gm::MODEL_MESH::CELL_POINT(int pos_cell)
{
	int flag = 0;
	std::vector<int> face_temp = cells[pos_cell].pos_faces;
	cells[pos_cell].pos_nodes = faces[face_temp[0]].pos_nodes;
	auto iter_end = cells[pos_cell].pos_nodes.cend() - 1;
	face_temp.erase(face_temp.cbegin());
	for (int j = 0; face_temp.size() != 1;)
	{
		if (*iter_end == faces[face_temp[j]].pos_nodes[0] || *iter_end == faces[face_temp[j]].pos_nodes[1])
		{
			if (*iter_end == faces[face_temp[j]].pos_nodes[0]) cells[pos_cell].pos_nodes.push_back(faces[face_temp[j]].pos_nodes[1]);
			else cells[pos_cell].pos_nodes.push_back(faces[face_temp[j]].pos_nodes[0]);
			iter_end = cells[pos_cell].pos_nodes.cend() - 1;
			face_temp.erase(face_temp.cbegin() + j);
			j = 0;
		}
		else j++;
	}
}

void gm::MODEL_MESH::DIRECT_POINTS(int pos_cell)
{
	NULL_VEC a, b;
	int pos_point_0 = cells[pos_cell].pos_nodes[0], pos_point_1 = cells[pos_cell].pos_nodes[1];
	a.x = points[pos_point_1].x - points[pos_point_0].x;
	a.y = points[pos_point_1].y - points[pos_point_0].y;
	b.x = cells[pos_cell].centre.x - points[pos_point_0].x;
	b.y = cells[pos_cell].centre.y - points[pos_point_0].y;
	if (DET(a, b) < 0)
	{
		for (int i = 0, j = (int)cells[pos_cell].pos_nodes.size() - 1; i != j && i < j; i++, j--)
		{
			int x = cells[pos_cell].pos_nodes[i];
			cells[pos_cell].pos_nodes[i] = cells[pos_cell].pos_nodes[j];
			cells[pos_cell].pos_nodes[j] = x;
		}
	}
}

void gm::MODEL_MESH::STRUCT_CELL()
{
	for (int pos_cell = 0; pos_cell < count_all_cells; pos_cell++)
	{
		CELL_POINT(pos_cell);
		CENTRE(cells[pos_cell].centre, cells[pos_cell].pos_nodes);
		RADIUS(cells[pos_cell].centre, cells[pos_cell].r);
		DIRECT_POINTS(pos_cell);
		SQUARE(pos_cell);
		CHORD(pos_cell);
	}
}

void gm::MODEL_MESH::STRUCT()
{
	STRUCT_FACE();
	STRUCT_CELL();
}

