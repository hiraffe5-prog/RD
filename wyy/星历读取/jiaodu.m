clear;clc
EARTH_emp_path = 'C:\Users\eric\Desktop\ToWYY\数据\WGC_FrameTransformation_20260312022956.csv';

% 读取欧拉角数据
[Shike_Line_Str, Euler_Angles] = nasa_csv_euler_read(EARTH_emp_path);

% 星历采样间隔
dt_ephem = 1;   % s

% 原始数据点数
Ttotal = size(Euler_Angles, 1);

% 原始时间轴
t_raw = 0:dt_ephem:(Ttotal-1)*dt_ephem;

% 欧拉角：度 -> 弧度
psi_raw   = Euler_Angles(:,1) * pi / 180;   % Angle 3
theta_raw = Euler_Angles(:,2) * pi / 180;   % Angle 2
phi_raw   = Euler_Angles(:,3) * pi / 180;   % Angle 1

% 三阶多项式拟合
psi_par   = polyfit(t_raw, psi_raw.', 3);
theta_par = polyfit(t_raw, theta_raw.', 3);
phi_par   = polyfit(t_raw, phi_raw.', 3);

% % 提取系数
% d_psi = psi_par(1);     c_psi = psi_par(2);     b_psi = psi_par(3);     a_psi = psi_par(4);
% d_theta = theta_par(1); c_theta = theta_par(2); b_theta = theta_par(3); a_theta = theta_par(4);
% d_phi = phi_par(1);     c_phi = phi_par(2);     b_phi = phi_par(3);     a_phi = phi_par(4);




