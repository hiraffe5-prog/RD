function [Shike_Line_Str,XYZ] = nasa_csv_read(EARTH_emp_path)
% Shike_Line_Str:N*-,string
% XYZ:N*3,double

%EARTH_emp_path = 'IAU_MARS_DE421_2021_12_10T04_00_00_000_2021_12_12T03_59_00_000.csv';
trace_E = readmatrix(EARTH_emp_path,'OutputType','string');
trace_E = trace_E(13:end-4,1:3);
Num_E = size(trace_E,1);
% 时间字符串单独存
Shike_Line_Str = strings(Num_E, 1);

% 坐标单独存，保持double精度
XYZ = zeros(Num_E, 3);

for ii = 1:Num_E
    % 按逗号拆分这一行
    t_3 = regexp(trace_E(ii,1), ',', 'split');

    % 提取X、Y、Z（第4、5、6个字段）
    XYZ(ii,1) = str2double(t_3(4));   % X (km)
    XYZ(ii,2) = str2double(t_3(5));   % Y (km)
    XYZ(ii,3) = str2double(t_3(6));   % Z (km)

    % 处理时间格式
    t_e_3 = regexp(t_3(1), '[TZ:. ]', 'split');
    t_mm = str2double(t_e_3(5));
    t_mm = sprintf('%03d', floor(t_mm/100));

    Shike_Line_Str(ii) = [char(t_e_3(1)), 'T', ...
                          char(t_e_3(2)), '-', ...
                          char(t_e_3(3)), '-', ...
                          char(t_e_3(4)), '-', ...
                          t_mm, 'Z'];
end
end