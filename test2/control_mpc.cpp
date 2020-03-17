#include"control_mpc.h"
#define PI 3.14159;
#define sign(d) (d < 0 ? -1 : d>0);
control_mpc::control_mpc(int* m_, int* p_, int* nwsr_, constraint_mpc* constraint_) :slover(*m_, *m_) {
	M = *m_;
	P = *p_;
	nWSR = *nwsr_;
	//初始化二次型
	result_u = MatrixXd::Zero(M, 1);
	res = result_u.data();//结果指针指向结果，两者之间可以相互变化
	H_qp = MatrixXd::Identity(M, M);
	g_qp = MatrixXd::Ones(M, 1);
	//生成约束矩阵
	lb_u = constraint_->min_u * MatrixXd::Ones(M, 1);
	ub_u = constraint_->max_u * MatrixXd::Ones(M, 1);
	A_cons = MatrixXd::Identity(M, M);
	A_cons.block(1, 0, M - 1, M - 1) -= MatrixXd::Identity(M - 1, M - 1);
	ubA_u = constraint_->max_deta_u * MatrixXd::Ones(M, 1);
	lbA_u = constraint_->min_deta_u * MatrixXd::Ones(M, 1);
	slover.init(H_qp.data(), g_qp.data(), A_cons.data(), lb_u.data(), ub_u.data(), lbA_u.data(), ubA_u.data(), nWSR);
	slover.getPrimalSolution(result_u.data());
	//初始化权值矩阵
	Q_qp = MatrixXd::Zero(4 * P, 4 * P);
	R_qp = MatrixXd::Zero(M - 1, M - 1);
	R_first = MatrixXd::Zero(1, 1);
	L_first = MatrixXd::Zero(1, M);
	L_first(0, 0) = 1;
	//确定Sx等数据大小
	Sx = MatrixXd::Zero(4 * P, 4);
	Su = MatrixXd::Zero(4 * P, M);
	Sc = MatrixXd::Zero(4 * P, 1);
}
void control_mpc::generate_leader(vector<float>&plan_x, vector<float>&plan_y) {
	leader= MatrixXd::Zero(plan_x.size(), 5);
	float pk[2], pk_1[2];
	//确定初始角度
	//vector<float>::iterator _x = pan_x->begin(), _y = pan_y->begin();
	pk[0] = plan_x[1] - plan_x[0], pk[1] = plan_y[1] - plan_y[0];
	leader(0, 3) = asin(pk[1] / pow(pk[0] * pk[0] + pk[1] * pk[1], 0.5));
	if (pk[0] < 0 && pk[1] < 0) {
		leader(0, 3) = -PI - leader(0, 3);
	}
	else if (pk[0] < 0 && pk[1]>0) {
		leader(0, 3) = PI - leader(0, 3);
	}
	else if (pk[0] < 0 && abs(pk[1]) < 0.0001) {
		leader(0, 3) = PI;
	}
	leader(0, 2) = 0, leader(0, 4) = 0;//第一个点s为0,kapa为0
	leader(0, 0) = plan_x[0], leader(0, 1) = plan_y[0];
	for (int i = 0; i < plan_x.size() - 2; i++) {
		pk[0] = plan_x[i+1] - plan_x[i], pk[1] = plan_y[i+1] - plan_y[i];
		pk_1[0] = plan_x[i + 2] - plan_x[i + 1], pk_1[1] = plan_y[i + 2] - plan_y[i + 1];
		leader(i + 1, 3) = leader(i, 3) + asin((pk[0] * pk_1[1] - pk_1[0] * pk[1]) / pow(pk[0] * pk[0] + pk[1] * pk[1], 0.5) / pow(pk_1[0] * pk_1[0] + pk_1[1] * pk_1[1], 0.5));
		leader(i + 1, 2) = leader(i, 2) + pow(pk[0] * pk[0] + pk[1] * pk[1], 0.5);
		leader(i + 1, 0) = plan_x[i + 1], leader(i + 1, 1) = plan_y[i + 1];
		leader(i + 1, 4) = (leader(i + 1, 3) - leader(i, 3)) / (leader(i + 1, 2) - leader(i, 2));
	}
	leader(plan_x.size() - 1, 0) = plan_x[plan_x.size() - 1], leader(plan_x.size() - 1, 1) = plan_y[plan_y.size() - 1];
	pk[0] = plan_x[plan_x.size() - 1] - plan_x[plan_x.size() - 2], pk[1] = plan_y[plan_y.size() - 1] - plan_y[plan_x.size() - 2];
	leader(plan_x.size() - 1, 2) = leader(plan_x.size() - 2, 2) + pow(pk[0] * pk[0] + pk[1] * pk[1], 0.5);
	leader(plan_x.size() - 1, 3) = leader(plan_x.size() - 2, 3);//最后一个点航向角不变；
	leader(plan_x.size() - 1, 4) = 0;
}
void control_mpc::sl_xy_change(car_state*_state,float* state_des_, double* leader_data, int leader_row) {//leader_data对应leader.data()
	double P_n[2], P_x_r[2], P_x_f[2];
	int g_cos_r, g_cos_f, g_sin_r;
	int s_g_later_f = 0, s_g_later_r = 0;
	double P_n_module;
	int i = 0;
	for (; i < leader_row - 1; i++) {
		//cout << i << endl;
		P_n[0] = leader_data[i + 1] - leader_data[i], P_n[1] = leader_data[leader_row + i - 1 + 1] - leader_data[leader_row + i - 1];
		P_x_r[0] = _state->x_real - leader_data[i], P_x_r[1] = _state->y_real - leader_data[leader_row + i - 1];
		g_cos_r = ((P_n[0] * P_x_r[0] + P_n[1] * P_x_r[1]));//判断是否为锐角，1为锐角-1为钝角，0垂直
		g_sin_r = sign(P_n[0] * P_x_r[1] - P_n[1] * P_x_r[0]);//判断是否重合
		P_x_f[0] = _state->x_real - leader_data[i + 1], P_x_f[1] = _state->y_real - leader_data[leader_row + i - 1 + 1];
		g_cos_f = sign(P_n[0] * P_x_f[0] + P_n[1] * P_x_f[1]);//判断下一个端点与端点向量之间的夹角是否为钝角。
		if ((g_sin_r == 0) & (((_state->x_real >= leader_data[i]) & (_state->x_real <= leader_data[i+1])))) {
			break;
		}
		else if ((g_cos_f * g_cos_r == -1) | (g_cos_r == 0)) {
			break;
		}
		//大圆角曲率
		else if (((s_g_later_f + s_g_later_r == 2) & (g_cos_f + g_cos_r == -2)) | ((s_g_later_f + s_g_later_r == -2) & (g_cos_f + g_cos_r == 2))) {
			break;
		}
		s_g_later_f = g_cos_f;
		s_g_later_r = g_cos_r;
	}
	P_n_module = sqrt(P_n[0] * P_n[0] + P_n[1] * P_n[1]);
	state_des_[0] = leader_data[i] + (P_n[0] * P_x_r[0] + P_n[1] * P_x_r[1]) / P_n_module;
	state_des_[1] = (P_n[0] * P_x_r[1] - P_n[1] * P_x_r[0]) / P_n_module;
	state_des_[2] = leader_data[2 * (leader_row - 1) + i] + (leader_data[2 * (leader_row - 1) + i + 1] - leader_data[2 * (leader_row - 1) + i]) * (state_des_[0] - leader_data[i]) / (leader_data[i + 1] - leader_data[i]);
	state_des_[3] = leader_data[3 * (leader_row - 1) + i] + (leader_data[3 * (leader_row - 1) + i + 1] - leader_data[3 * (leader_row - 1) + i]) * (state_des_[0] - leader_data[i]) / (leader_data[i + 1] - leader_data[i]);
}
void control_mpc::generate_weight_mat(weight_mpc* weight_) {
	for (int i = 0; i < P; i++) {
		Q_qp(4*i, 4*i) = weight_->weight_eg;
		Q_qp(4 * i + 1, 4 * i + 1) = weight_->weight_eg_dot;
		Q_qp(4 * i + 2, 4 * i + 2) = weight_->weight_hg;
		Q_qp(4 * i + 3, 4 * i + 3) = weight_->weight_hg_dot;
	}
	for (int i = 0; i < M - 1; i++) {
		R_qp(i, i) = weight_->weight_detaU;
	}
	R_first(0, 0) = weight_->weight_detaU_first;
}
void control_mpc::update_state_martix(car_state*_state) {
	float state_des[4];//期望状态[s_des,l_des,heading_des,kapa_des]
	sl_xy_change(_state, state_des, leader.data(), leader.size() / 5);
	Matrix<double, 4, 1> state_erro;
	state_erro(0, 0) = state_des[1], state_erro(1, 0) = _state->vy + _state->vx * sin(_state->heading - state_des[2]);
	state_erro(2, 0) = _state->heading - state_des[2], state_erro(3, 0) = _state->heading_rate - _state->vx * state_des[3];
	Ac(0, 1) = 1;
	Ac(1, 1) = -(Cr + Cf) / mass / _state->vx;
	Ac(1, 2) = (Cr + Cf) / mass;
	Ac(1, 3) = (Lr * Cr - Lf * Cf) / mass / _state->vx;
	Ac(2, 3) = 1;
	Ac(3, 1) = (Lr * Cr - Lf * Cf) / Iz / _state->vx;
	Ac(3, 2) = (Lf * Cf - Lr * Cr) / Iz; 
	Ac(3, 3) = -(Lf * Lf * Cf + Lr * Lr *Cr) / Iz / _state->vx;
	Bc(1, 0) = Cf / mass, Bc(3, 0) = Lf * Cf / Iz;
	Cc(1, 0) = (Lr * Cr - Lf * Cf) / mass * state_des[3] - (_state->vx) *(_state->vx) * state_des[3];
	Cc(3, 0) = -(Lf * Lf * Cf - Lr * Lr * Cr) / Iz * state_des[3];
	Matrix<double, 4, 4>Ac_inv = (MatrixXd::Identity(4, 4) - Ts / 2 * Ac).inverse();
	Ac =  Ac_inv* (MatrixXd::Identity(4, 4) + Ts / 2 * Ac);
	Bc = Ts * Ac_inv * Bc;
	Cc = Ts * Ac_inv * Cc;
	Matrix<double, 4, 4> Sx_A;
	Matrix<double, 4, 1> Sc_AC;
	Sx_A = Ac; Sc_AC = Cc;
	Sx.middleRows(0, 4) = Sx_A;
	Su.block(0, 0, 4, 1) = Bc;
	Sc.middleRows(0, 4) = Sc_AC;
	for (int i = 1; i < P; ++i) {
		Su.block(i * 4, 0, 4, 1) = Sx_A * Bc;
		for (int j = 1; j <= i; ++j) {
			Su.block(i * 4, j, 4, 1) = Su.block(i * 4 - 4, j - 1, 4, 1);
		}
		Sc_AC += Sx_A * Cc;
		Sc.middleRows(i * 4, 4) = Sc_AC;
		Sc.middleRows(i * 4, 4) = Sx_A * Cc;
		Sx_A *= Ac;
		Sx.middleRows(4 * i, 4) = Sx_A;
	}
	cout << Sx << endl << Su << endl << Sc << endl;
	H_qp = 2 * (Su.transpose() * Q_qp * Su + A_cons.transpose() * R_qp * A_cons + L_first.transpose() * R_first * L_first);
	g_qp = 2 * (Sx.transpose() * Q_qp * state_erro + Su.transpose() * Q_qp * Su - L_first.transpose() * R_first * _state->u_k_1);
	ubA_u(0, 0) += _state->u_k_1;
	lbA_u(0, 0) += _state->u_k_1;
}
void control_mpc::control_mpc_slover() {
	slover.hotstart(H_qp.data(), g_qp.data(), A_cons.data(), lb_u.data(), ub_u.data(), lbA_u.data(), ubA_u.data(), nWSR);
	slover.getPrimalSolution(result_u.data());
}