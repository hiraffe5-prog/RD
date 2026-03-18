clc; clear;

%% 1. 定义符号变量
syms X Y Z t2 real
syms a b c d e f R real

% 三个欧拉角三次多项式系数
syms a_psi b_psi c_psi d_psi real
syms a_theta b_theta c_theta d_theta real
syms a_phi b_phi c_phi d_phi real

%% 2. 定义欧拉角关于 t2 的表达式
psi   = a_psi   + b_psi*t2   + c_psi*t2^2   + d_psi*t2^3;
theta = a_theta + b_theta*t2 + c_theta*t2^2 + d_theta*t2^3;
phi   = a_phi   + b_phi*t2   + c_phi*t2^2   + d_phi*t2^3;

%% 3. 定义三个旋转矩阵
Rz_psi = [ cos(psi),  sin(psi), 0;
          -sin(psi),  cos(psi), 0;
                0,         0,   1];

Rx_theta = [1,      0,           0;
            0, cos(theta),  sin(theta);
            0,-sin(theta),  cos(theta)];

Rz_phi = [ cos(phi),  sin(phi), 0;
          -sin(phi),  cos(phi), 0;
               0,         0,    1];

%% 4. 定义转换矩阵 Tran(t2)
Tran = Rz_psi * Rx_theta * Rz_phi;

%% 5. 定义目标点坐标 pt(t2)
pt = Tran * [X; Y; Z];

%% 6. 定义雷达坐标 pr(t1), pr(t3)
pr1 = [a; b; c];
pr3 = [d; e; f];

%% 7. 构造 F1
F1 = sqrt((pr1(1)-pt(1))^2 + (pr1(2)-pt(2))^2 + (pr1(3)-pt(3))^2) ...
   + sqrt((pt(1)-pr3(1))^2 + (pt(2)-pr3(2))^2 + (pt(3)-pr3(3))^2) ...
   - R;
%% 8. 求 F1 对 X 的偏导数
dF1_dX = diff(F1, X);

%% 9. 化简结果
dF1_dX = simplify(dF1_dX);

%% 10. 显示结果
disp('F1 对 X 的偏导数为：');
disp(dF1_dX)
