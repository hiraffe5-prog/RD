clear;clc
EARTH_emp_path = 'C:\Users\eric\Desktop\ToWYY\数据\WGC_StateVector_20260312022150.csv';

% 读取雷达在 J2000 月心惯性系（MCI）下的星历
[Shike_Line_Str, XYZ] = nasa_csv_read(EARTH_emp_path);

% 三维坐标（km -> m）
x_Radar_MCI = XYZ(:,1) * 1000;
y_Radar_MCI = XYZ(:,2) * 1000;
z_Radar_MCI = XYZ(:,3) * 1000;

% 星历采样间隔
dt_ephem = 1;              % s

% 雷达脉冲重复时间
PRT = 10e-3;               % s

% 原始星历点数
Ttotal = length(x_Radar_MCI);

% 原始星历时间轴（单位：s）
t_raw = 0:dt_ephem:(Ttotal-1)*dt_ephem;

% 对三个坐标分量分别进行三阶多项式拟合
x_par_radar = polyfit(t_raw, x_Radar_MCI, 3);   % [d_x, c_x, b_x, a_x]
y_par_radar = polyfit(t_raw, y_Radar_MCI, 3);   % [d_y, c_y, b_y, a_y]
z_par_radar = polyfit(t_raw, z_Radar_MCI, 3);   % [d_z, c_z, b_z, a_z]

% 总观测时长
T_end = t_raw(end);

% 雷达慢时间轴
ta = 0:PRT:T_end;

% 在慢时间轴上计算拟合后的雷达轨迹
x_radar_ta = polyval(x_par_radar, ta);
y_radar_ta = polyval(y_par_radar, ta);
z_radar_ta = polyval(z_par_radar, ta);

% 组合成雷达在慢时间轴上的三维位置
p_r_ta = [x_radar_ta(:), y_radar_ta(:), z_radar_ta(:)];




