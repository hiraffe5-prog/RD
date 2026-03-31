% BP对月成像处理
clc
if isempty(gcp('nocreate'))
    parpool;
end
addpath 'C:\Users\eric\Desktop\月球面目标仿真代码\子函数'
addpath 'C:\Users\eric\Desktop\wyy_RD定位\星历读取'
%% 参数设置
c = 299792458;
R = 1737400;          % 月球参考半径(m)
load('C:\Users\eric\Desktop\月球面目标仿真代码\图像数据\amsum_new.mat')
Tcoh = 120;
PRF  = 20;
fc   = 3.1e9;
Tp   = 1e-4;          % 脉冲持续时间
B    = 0.5e6;
k    = B / Tp;
fs   = 2 * B;
Tr   = 1 / PRF;
Na   = round(PRF * Tcoh);
Na   = Na + mod(Na, 2);   % 保证Na为偶数
Tcoh = Na * Tr;
lambda = c / fc;

% lon  = -11;
% lat  = -43;
% 
% center_x = (R-4470).*cosd(lat).*cosd(lon);  %散射体的x坐标
% center_y = (R-4470).*cosd(lat).*sind(lon);  %散射体的y坐标
% center_z = (R-4470).*sind(lat); 

%% 坐标拟合
EARTH_emp_path_MCI = 'C:\Users\eric\Desktop\月球面目标仿真代码\输入文件\WGC_StateVector_20260327093305_radar_J2000.csv';

% 读取雷达在 J2000 月心惯性系（MCI）下的星历
[~, XYZ_MCI] = nasa_csv_read(EARTH_emp_path_MCI);

% 三维坐标（km -> m）
x_Radar_MCI = XYZ_MCI(:,1) * 1000;
y_Radar_MCI = XYZ_MCI(:,2) * 1000;
z_Radar_MCI = XYZ_MCI(:,3) * 1000;

% 星历采样间隔
dt_ephem = 1;              % s

% 原始星历点数
Ttotal = length(x_Radar_MCI);

% 原始星历时间轴（单位：s）
t_raw = 0:dt_ephem:(Ttotal-1)*dt_ephem;

% 对三个坐标分量分别进行三阶多项式拟合
x_par_radar = polyfit(t_raw, x_Radar_MCI, 3);   % [d_x, c_x, b_x, a_x]
y_par_radar = polyfit(t_raw, y_Radar_MCI, 3);   % [d_y, c_y, b_y, a_y]
z_par_radar = polyfit(t_raw, z_Radar_MCI, 3);   % [d_z, c_z, b_z, a_z]

% 雷达慢时间轴
ta = (0:Na-1)*Tr;

% 在慢时间轴上计算拟合后的雷达轨迹
x_radar_ta = polyval(x_par_radar, ta);
y_radar_ta = polyval(y_par_radar, ta);
z_radar_ta = polyval(z_par_radar, ta);

% 组合成雷达在慢时间轴上的三维位置
p_r_ta = [x_radar_ta(:), y_radar_ta(:), z_radar_ta(:)];

%% 欧拉角拟合
Euler_emp_path = 'C:\Users\eric\Desktop\月球面目标仿真代码\输入文件\WGC_FrameTransformation_20260327104919.csv';

% 读取欧拉角数据
[~, Euler_Angles] = nasa_csv_euler_read(Euler_emp_path);

% 欧拉角：度 -> 弧度
psi_raw   = Euler_Angles(:,1) * pi / 180;   % Angle 3
theta_raw = Euler_Angles(:,2) * pi / 180;   % Angle 2
phi_raw   = Euler_Angles(:,3) * pi / 180;   % Angle 1

% 三阶多项式拟合
psi_par   = polyfit(t_raw, psi_raw.', 3);
theta_par = polyfit(t_raw, theta_raw.', 3);
phi_par   = polyfit(t_raw, phi_raw.', 3);

% 提取系数
d_psi = psi_par(1);     c_psi = psi_par(2);     b_psi = psi_par(3);     a_psi = psi_par(4);
d_theta = theta_par(1); c_theta = theta_par(2); b_theta = theta_par(3); a_theta = theta_par(4);
d_phi = phi_par(1);     c_phi = phi_par(2);     b_phi = phi_par(3);     a_phi = phi_par(4);

lon_range=[-13.5,-8.5];
lat_range=[-40.5,-45.5];

lon_lat_range=[lon_range,lat_range];
del_lat = 0.005;
del_lon = 0.005;
lon = linspace(lon_lat_range(1), lon_lat_range(2), abs(lon_lat_range(2)-lon_lat_range(1))/del_lon);
lat = linspace(lon_lat_range(3), lon_lat_range(4), abs(lon_lat_range(4)-lon_lat_range(3))/del_lat);

globalDEM_path = 'C:\Users\eric\Desktop\wyy_RD定位\读取DEM\DEM_global_moon.mat';
[DEM_1]=DEM_crop_fun(globalDEM_path,lon_lat_range); 
[DEM_interp]=DEM_resample(DEM_1,lon_lat_range,lon,lat);

lon_idx = (lon >= -13) & (lon <= -9);
lat_idx = (lat >= -45) & (lat <= -41);
lon_crop = lon(lon_idx);
lat_crop = lat(lat_idx);

amsum_crop = amsum_new(lat_idx, lon_idx);

% 裁剪DEM
DEM_crop = DEM_interp(lat_idx, lon_idx);

% figure;
% imagesc(lon_crop, lat_crop, DEM_crop)
% set(gca,'YDir','normal')
% xlabel('经度/deg')
% ylabel('纬度/deg')
% title('裁剪后的DEM')

[R_axis, fd_axis, Rm, fd, T2m, TarXYZm] = BP_to_RD_nonstopgo( ...
    amsum_crop, DEM_crop, lon_crop, lat_crop, Na, Tr, ...
    x_par_radar, y_par_radar, z_par_radar, ...
    a_psi, b_psi, c_psi, d_psi, ...
    a_theta, b_theta, c_theta, d_theta, ...
    a_phi, b_phi, c_phi, d_phi, ...
    lambda, R, c);

    % 绘制反投到RD域的bp成像结果
    % figure; scatter(Rm(:),Vm(:),[],abs(amsum_new(:)),'.');
    figure;
    scatter(Rm(:), fd(:), 1, abs(amsum_crop(:)), 'filled');
    xlabel('双程斜距 R (m)');
    ylabel('多普勒 f_d (Hz)');
    title('BP结果映射到RD域');
    % colormap jet;
    grid on;
    box on;
    % 设置坐标范围
    xlim([min(Rm(:)), max(Rm(:))]);
    ylim([min(fd(:)), max(fd(:))]);
