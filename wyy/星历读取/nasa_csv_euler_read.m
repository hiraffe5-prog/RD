function [Shike_Line_Str,Euler_Angles] = nasa_csv_euler_read(EARTH_emp_path)
% Shike_Line_Str : N*1 string
% Euler_Angles   : N*3 double
%                  第1列 Angle 3 = psi
%                  第2列 Angle 2 = theta
%                  第3列 Angle 1 = phi

%EARTH_emp_path = 'IAU_MARS_DE421_2021_12_10T04_00_00_000_2021_12_12T03_59_00_000.csv';
trace_E = readmatrix(EARTH_emp_path,'OutputType','string');
trace_E = trace_E(1:end-4,1:4);

Num_E = size(trace_E,1);
% 角度直接整体转换
Euler_Angles = str2double(trace_E(:,2:4));
% 时间单独转换
Shike_Line_Str = strings(Num_E,1);
for ii = 1:Num_E
    t_e_3 = regexp(trace_E(ii,1),'[TZ:. ]', 'split');
    t_mm = str2double(t_e_3(5));
    t_mm = sprintf('%03d',floor(t_mm/100));
    Shike_Line_Str(ii) = [char(t_e_3(1)),'T',char(t_e_3(2)),'-', ...
                          char(t_e_3(3)),'-',char(t_e_3(4)),'-',t_mm,'Z'];
end

end