clc; clear;

%% 1. 定义符号变量
syms X Y Z t t2 real
syms a b c d e f real          % pr(t1), pr(t3)
syms g h i j k l real          % pr_dot(t1), pr_dot(t3)
syms lambda fd real

% 欧拉角三次多项式系数
syms a_psi b_psi c_psi d_psi real
syms a_theta b_theta c_theta d_theta real
syms a_phi b_phi c_phi d_phi real

%% 2. 定义欧拉角关于一般时间 t 的表达式
psi_t   = a_psi   + b_psi*t   + c_psi*t^2   + d_psi*t^3;
theta_t = a_theta + b_theta*t + c_theta*t^2 + d_theta*t^3;
phi_t   = a_phi   + b_phi*t   + c_phi*t^2   + d_phi*t^3;

%% 3. 定义旋转矩阵 Tran(t)
Rz_psi = [ cos(psi_t),  sin(psi_t), 0;
          -sin(psi_t),  cos(psi_t), 0;
                 0,          0,     1];

Rx_theta = [1,      0,           0;
            0, cos(theta_t),  sin(theta_t);
            0,-sin(theta_t),  cos(theta_t)];

Rz_phi = [ cos(phi_t),  sin(phi_t), 0;
          -sin(phi_t),  cos(phi_t), 0;
               0,           0,      1];

Tran_t = Rz_psi * Rx_theta * Rz_phi;

%% 4. 构造 pt(t)
pt_t = Tran_t * [X; Y; Z];

%% 5. 对 t 求导，得到 pt_dot(t)
dpt_t = diff(pt_t, t);

%% 6. 在 t = t2 处取值
pt_t2  = subs(pt_t,  t, t2);
dpt_t2 = subs(dpt_t, t, t2);

%% 7. 已知量向量
pr1  = [a; b; c];
pr3  = [d; e; f];
vpr1 = [g; h; i];
vpr3 = [j; k; l];

%% 8. 构造 F2
term1 = ( (pt_t2 - pr1).' * (vpr1 - dpt_t2) ) / (lambda * sqrt((pt_t2(1)-a)^2 + (pt_t2(2)-b)^2 + (pt_t2(3)-c)^2));
term2 = ( (pr3 - pt_t2).' * (dpt_t2 - vpr3) ) / (lambda * sqrt((pt_t2(1)-d)^2 + (pt_t2(2)-e)^2 + (pt_t2(3)-f)^2));

F2 = simplify(term1 + term2 - fd);

disp('F2 = ')
disp(F2)