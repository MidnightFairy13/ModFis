#pragma once
#include "Model.h"
#include "Values.h"
namespace gm
{
	class FISICS
	{
	public:
		MODEL_MESH Mesh;
		Inf_Cond inf;
		std::vector<OPT_CELL> Cells;
		std::vector<OPT_FACE> Faces;
		FISICS();
		void CALCULATING();
		std::string INF_FLOW(int pos_flow, bool fluid = 1);
		double Tet();
		double dGi(int pos_zone, int pos_flow, bool fluid = 1);
		std::vector<double> dGi(int pos_zone, bool fluid = 1);
	private:
		double Fi(int pos_face, int pos_flow, bool fluid = 1);
		std::vector<double> Fi(int pos_face, bool fluid = 1);
		void START_FLUID(gm::OPT_CELL& cell);
		void START_PART(gm::OPT_CELL& cell);
		void MAX_DIF(int pos_cell, double& U, double& V);
		double tau();
		double Un(int pos_cell, int pos_face, bool fluid = 1);
		double Ut(int pos_cell, int pos_face, bool fluid = 1);
		void T_1(double& f1, double& f2, int pos_face);
		void T_1(int pos_face, bool fluid = 1);
		void F_FLUID(int pos_face, double _U, double _V, double _P, double _Ro, bool inlet = 1);
		std::vector<double> F_FLUID(int pos_cell, int pos_face);
		void F_PART(int pos_face, double _U, double _V, double _T, double _Ro);
		std::vector<double> F_PART(int pos_cell, int pos_face);
		double dG(int pos_flow, int iter, bool fluid = 1);
		std::vector<double> dG(int iter, bool fluid = 1);

		std::string FLAG_FLOW(int iter, std::vector<double>& max, bool fluid = 1);
		bool WORK_STOP(int count, int& finish, std::string stop_string);
		bool STOP();
		bool FINISH(int count, int& finish);
		bool FLOW_STOP(std::string flow_stop);
		void WRITE();

		void FLOW_FACE();
		void FLOW_FACE(int pos_zone);
		void FLOW_PRESSURE_OUTLET(int pos_zone);//thread
		void FLOW_MASS_FLOW_INLET(int pos_zone);//thread
		void FLOW_VELOCITY_INLET(int pos_zone);//thread
		void FLOW_OUTFLOW(int pos_zone);//thread
		void FLOW_VAN_LEER(int pos_zone);//thread

		void PRESSURE_OUTLET(int pos_face);
		double MFI(int pos_face, double Pi);
		void MASS_FLOW_INLET(int pos_face);
		void VELOCITY_INLET(int pos_face);
		void OUTFLOW(int pos_face);
		bool WALL(int pos_face);
		void VAN_LEER(int pos_face);
		std::vector<double> VAN_LEER(int pos_cell, int pos_face, int dir = 1);
		void FLOW_PART(int pos_face);

		void FLOW_CELL(double dt);//thread
		void FLOW_CELL_MAIN(double dt, int pos_cell);
		void FLOW_CELL(double dt, int pos_cell, bool fluid = 1);
		std::vector<double> FLOW_CELL(int pos_cell, bool fluid = 1);
		std::vector<double> AXIS(int pos_cell, bool fluid = 1);
		std::vector<double> AXIS_FLUID(int pos_cell);
		std::vector<double> AXIS_PART(int pos_cell);
		std::vector<double> BURN(int pos_cell, bool fluid = 1);
		void INTERACTION(double dt, int pos_cell);

		void Yi(int pos_cell, bool fluid = 1);
		void NEW_OPT(int pos_cell, bool fluid = 1);
		void NEW_OPT_UV(int pos_cell, bool fluid = 1);
	};
}