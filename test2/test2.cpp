#include"control_mpc.h"
#include <iostream>
using namespace std;
int main() {
	float _x_real = 0, _y_real = 0, _heading = 0, _u_k_1 = 0, _vx = 50, _vy = 0, _heading_rate = 0;//获得实时车辆状态
	car_state _state_current(&_x_real, &_y_real, &_heading, &_u_k_1, &_vx, &_vy, &_heading_rate);
	weight_mpc _weight;//设置权重
	constraint_mpc _constraint(1.0f, 1.0f, 1.0f, 1.0f);//设置约束
	vector<float>_plan_x, _plan_y;//输入规划数据,实时更新
	//需更改private数据
	int _m = 4, _p = 8, _nwsr = 20;
	control_mpc _Control(&_m, &_p, &_nwsr, &_constraint);
	_Control.generate_weight_mat(&_weight);
	cout << _plan_x.size() << endl;
	while (_plan_x.size() != 0 && _plan_y.size() != 0) {
		_Control.generate_leader(_plan_x, _plan_y);
		_Control.update_state_martix(&_state_current);
		_Control.control_mpc_slover();
		cout << *_Control.res << endl;
	}
	return 0;
}