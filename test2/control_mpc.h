#include <iostream>
#include "qpOASES.hpp"
#include<Eigen/Dense>
#include<vector>
using namespace Eigen;
//using namespace qpOASES;
using namespace std;
struct car_state {
	float x_real, y_real;
	float heading, u_k_1, heading_rate;
	float vx, vy;
	car_state(float* x_, float* y_, float* h_, float* u_,float*vx_,float*vy_,float*h_rate) :x_real(*x_), y_real(*y_), heading(*h_), u_k_1(*u_),vx(*vx_),vy(*vy_),heading_rate(*h_rate) {};
};
struct weight_mpc {
	float weight_eg, weight_eg_dot, weight_hg, weight_hg_dot;
	float weight_detaU, weight_detaU_first;
	weight_mpc(float _eg = 1.0f,float _eg_dot = 1.0f,float _hg = 1.0f, float _hg_dot = 1.0f, float _detaU = 1.0f, float _detaU_first = 1.0f) {
		weight_eg = _eg;
		weight_eg_dot = _eg_dot;
		weight_hg = _hg;
		weight_hg_dot = _hg_dot;
		weight_detaU = _detaU;
		weight_detaU_first = _detaU_first;
	}
};
struct constraint_mpc {
	float max_deta_u, min_deta_u;
	float max_u, min_u;
	constraint_mpc(float max_deta_,float min_deta_,float max_,float min_):max_deta_u(max_deta_),min_deta_u(min_deta_),max_u(max_),min_u(min_){}
};
class control_mpc {
public:
	control_mpc(int* m_, int* p_, int* nwsr_,  constraint_mpc* constraint_);
	void generate_weight_mat(weight_mpc* weight_);
	void update_state_martix(car_state* _state);
	Matrix<double, Dynamic, 5> leader;//[x,y,s,heading,kapa]
	void generate_leader(vector<float>& plan_x, vector<float>& plan_y);
	void sl_xy_change(car_state* _state, float* state_des_, double* leader_data, int leader_row);
	void control_mpc_slover();
	double* res;//Ö¸Ïòresult_u;
private:
	qpOASES::SQProblem slover;
	int M, P, nWSR;
	Matrix<double,4,4> Ac;
	Matrix<double, 4, 1> Bc, Cc, state;
	Matrix<double, Dynamic, 4> Sx;
	Matrix<double, Dynamic, Dynamic> Su;
	Matrix<double, Dynamic, 1> Sc;
	Matrix<double, Dynamic, Dynamic> H_qp;
	Matrix<double, Dynamic, 1> g_qp;
	Matrix<double, Dynamic, Dynamic> Q_qp, R_qp, R_first;
	Matrix<double, Dynamic, 1> lb_u, ub_u;
	Matrix<double, Dynamic, Dynamic> A_cons, L_first;
	Matrix<double, Dynamic, 1> lbA_u, ubA_u;
	Matrix<double, Dynamic, 1> result_u;
	double Cf = 1, Cr = 1, mass = 1, Iz = 1, Lr = 1, Lf = 1;
	float Ts = 0.7;
};