clear;clc;
% 研究目标的经纬度（为了简化计算过程，真正的场景中心会被设置为截取的DEM数据矩阵的中心，并不会完全等于这个取值）
target = [-11,-43];                                                 % 研究目标的经纬度【经度，纬度】（°）
globalDEM_path = './DEM_global_moon.mat';
%% 研究目标位置附近DEM数据截取
% 截取经纬度范围
lon_range = [target(1)-2,target(1)+2];      % 截取经度范围（°）
lat_range = [target(2)+2,target(2)-2];      % 截取纬度范围（°）
% 构建经纬度坐标轴
lon_lat_range=[lon_range,lat_range];

% del_lat = 0.05;
% del_lon = 0.05;
% lon = linspace(lon_lat_range(1), lon_lat_range(2), round(abs(lon_lat_range(2)-lon_lat_range(1))/del_lon)+1);
% lat = linspace(lon_lat_range(3), lon_lat_range(4), round(abs(lon_lat_range(4)-lon_lat_range(3))/del_lat)+1);
% % 构建二维经纬度网格
% [lon_2D,lat_2D] = meshgrid(lon,lat);

% 截取DEM数据
DEM = DEM_crop_fun(globalDEM_path,lon_lat_range);
% 读取DEM中心点经纬度和高程
[center_lon, center_lat, center_h, row_c, col_c, lon1, lat1] = get_dem_center(DEM, lon_lat_range);

% % 获取任意点的高程
lon_q=-10;lat_q=-43;
h_interp = interp2(lon1, lat1, DEM, lon_q, lat_q, 'linear');


