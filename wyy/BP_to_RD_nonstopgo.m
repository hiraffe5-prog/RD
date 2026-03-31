function [R_axis, fd_axis, Rm, fd, T2m, TarXYZm] = BP_to_RD_nonstopgo( ...
    y_bp, DEM_interp, lon, lat, Na, Tr, ...
    x_par_radar, y_par_radar, z_par_radar, ...
    a_psi, b_psi, c_psi, d_psi, ...
    a_theta, b_theta, c_theta, d_theta, ...
    a_phi, b_phi, c_phi, d_phi, ...
    lambda, Rmoon, c0, dR_out, dfd_out)
% ============================================================
% 功能：
%   将 BP 经-纬度网格复图像 y_bp 投影到非停走 RD 域
%
% 输入：
%   y_bp        : BP复图像，大小 Ny x Nx
%   DEM_interp  : 与 y_bp 同网格的高程矩阵，Ny x Nx
%   lon, lat    : 经度/纬度向量（单位：deg）
%   Na, Tr      : 脉冲数、PRT
%   x/y/z_par_radar : 雷达MCI坐标三次拟合系数 [1x4]
%   a_psi...d_phi   : 欧拉角三次拟合系数
%   lambda      : 波长
%   Rmoon       : 月球参考半径（m）
%   c0          : 光速
%   dR_out      : RD输出距离轴间隔（m），可传 [] 自动设置
%   dfd_out     : RD输出多普勒轴间隔（Hz），可传 [] 自动设置
%
% 输出：
%   y_rd        : RD域复图像
%   R_axis      : 距离轴（双程斜距，m）
%   fd_axis     : 多普勒轴（Hz）
%   Rm          : 每个BP像素对应的非停走双程斜距，Ny x Nx
%   Vm          : 每个BP像素对应的非停走多普勒，Ny x Nx
%
% 说明：
%   1) 参考时刻取中间慢时间：t1 = (round(Na/2)+1 - 1)*Tr
%   2) Rm 保存为双程斜距，不除以2
%   3) 多普勒采用你给出的非停走公式对应写法
%   4) RD投影时保留复数，不取 abs(y_bp)
% ============================================================

if nargin < 24 || isempty(c0)
    c0 = 299792458;
end
if nargin < 23 || isempty(Rmoon)
    Rmoon = 1737400;
end

[Ny, Nx] = size(y_bp);
assert(all(size(DEM_interp) == [Ny, Nx]), 'DEM_interp 尺寸必须与 y_bp 一致');

% 中间慢时间索引
midIdx = round(Na/2) + 1;
t1 = (midIdx - 1) * Tr;

% 经纬网格
[lon2D, lat2D] = meshgrid(lon, lat);

% 经纬高 -> MCMF
X_mcmf = (Rmoon + DEM_interp) .* cosd(lat2D) .* cosd(lon2D);
Y_mcmf = (Rmoon + DEM_interp) .* cosd(lat2D) .* sind(lon2D);
Z_mcmf = (Rmoon + DEM_interp) .* sind(lat2D);

%% ---------------- 逐像素计算双程斜距与多普勒 ----------------
Rm = zeros(Ny, Nx);   % 双程斜距
fd = zeros(Ny, Nx);   % 非停走多普勒
T2m = zeros(Ny, Nx);   % 每个像素对应的 t2
TarXYZm = zeros(Ny, Nx, 3);   % 每个像素的MCMF坐标


fprintf('开始计算每个BP像素的非停走双程斜距和多普勒...\n');

parfor iy = 1:Ny
    Rm_row = zeros(1, Nx);
    fd_row = zeros(1, Nx);
    T2_row = zeros(1, Nx);
    Tar_row  = zeros(Nx, 3); 

    for ix = 1:Nx
        Tar_xyz = [X_mcmf(iy,ix); Y_mcmf(iy,ix); Z_mcmf(iy,ix)];
        Tar_row(ix, :) = Tar_xyz(:).';
        [R_two_way, fd_val, t2_val] = calc_one_pixel_R_fd_nonstop( ...
            Tar_xyz, t1, ...
            x_par_radar, y_par_radar, z_par_radar, ...
            a_psi, b_psi, c_psi, d_psi, ...
            a_theta, b_theta, c_theta, d_theta, ...
            a_phi, b_phi, c_phi, d_phi, ...
            lambda, c0);

        Rm_row(ix) = R_two_way;
        fd_row(ix) = fd_val;
        T2_row(ix) = t2_val;

    end

    Rm(iy,:) = Rm_row;
    fd(iy,:) = fd_row;
    T2m(iy,:)= T2_row;
    TarXYZm(iy,:,:) = Tar_row;
    
    if mod(iy, max(1, floor(Ny/50))) == 0 || iy == Ny
        fprintf('RD几何计算进度：%d / %d (%.2f%%)\n', iy, Ny, 100*iy/Ny);
    end
end
end

%% ============================================================
% 单像素：非停走双程斜距 + 非停走多普勒
%% ============================================================
function [R_two_way, fd_val, t2] = calc_one_pixel_R_fd_nonstop( ...
    Tar_xyz, t1, ...
    x_par_radar, y_par_radar, z_par_radar, ...
    a_psi, b_psi, c_psi, d_psi, ...
    a_theta, b_theta, c_theta, d_theta, ...
    a_phi, b_phi, c_phi, d_phi, ...
    lambda, c0)

Tar_xyz = Tar_xyz(:);

%% 欧拉角
psi_fun   = @(t) d_psi   * t.^3 + c_psi   * t.^2 + b_psi   * t + a_psi;
theta_fun = @(t) d_theta * t.^3 + c_theta * t.^2 + b_theta * t + a_theta;
phi_fun   = @(t) d_phi   * t.^3 + c_phi   * t.^2 + b_phi   * t + a_phi;

dpsi_fun   = @(t) 3*d_psi   * t.^2 + 2*c_psi   * t + b_psi;
dtheta_fun = @(t) 3*d_theta * t.^2 + 2*c_theta * t + b_theta;
dphi_fun   = @(t) 3*d_phi   * t.^2 + 2*c_phi   * t + b_phi;

Rz = @(ang) [ cos(ang),  sin(ang), 0;
             -sin(ang),  cos(ang), 0;
                   0   ,      0   , 1 ];

Rx = @(ang) [1,     0     ,      0    ;
             0, cos(ang),  sin(ang);
             0, -sin(ang), cos(ang)];

dRz = @(ang,dang) [ -sin(ang)*dang,  cos(ang)*dang, 0;
                    -cos(ang)*dang, -sin(ang)*dang, 0;
                          0      ,        0       , 0 ];

dRx = @(ang,dang) [0, 0, 0;
                   0, -sin(ang)*dang,  cos(ang)*dang;
                   0, -cos(ang)*dang, -sin(ang)*dang];

% MCMF -> MCI
tran = @(t) Rz(psi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t));

% 时间导数
dtran = @(t) ...
    dRz(psi_fun(t), dpsi_fun(t)) * Rx(theta_fun(t)) * Rz(phi_fun(t)) + ...
    Rz(psi_fun(t)) * dRx(theta_fun(t), dtheta_fun(t)) * Rz(phi_fun(t)) + ...
    Rz(psi_fun(t)) * Rx(theta_fun(t)) * dRz(phi_fun(t), dphi_fun(t));

%% 雷达MCI位置与速度
xr_fun = @(t) x_par_radar(1)*t.^3 + x_par_radar(2)*t.^2 + x_par_radar(3)*t + x_par_radar(4);
yr_fun = @(t) y_par_radar(1)*t.^3 + y_par_radar(2)*t.^2 + y_par_radar(3)*t + y_par_radar(4);
zr_fun = @(t) z_par_radar(1)*t.^3 + z_par_radar(2)*t.^2 + z_par_radar(3)*t + z_par_radar(4);

vxr_fun = @(t) 3*x_par_radar(1)*t.^2 + 2*x_par_radar(2)*t + x_par_radar(3);
vyr_fun = @(t) 3*y_par_radar(1)*t.^2 + 2*y_par_radar(2)*t + y_par_radar(3);
vzr_fun = @(t) 3*z_par_radar(1)*t.^2 + 2*z_par_radar(2)*t + z_par_radar(3);

pr_fun  = @(t) [xr_fun(t);  yr_fun(t);  zr_fun(t)];
vpr_fun = @(t) [vxr_fun(t); vyr_fun(t); vzr_fun(t)];

%% 目标MCI位置与速度
pt_fun  = @(t) tran(t)  * Tar_xyz;
vpt_fun = @(t) dtran(t) * Tar_xyz;

%% 已知 t1 的雷达位置
pr_t1  = pr_fun(t1);
vpr_t1 = vpr_fun(t1);

%% ---------- 求 t2 ----------
fun_t2 = @(t) norm(pt_fun(t) - pr_t1) - c0 * (t - t1);

R1_init = norm(pt_fun(t1) - pr_t1);
t2_init = t1 + R1_init / c0;

t2 = local_fzero_with_expand(fun_t2, t2_init, t1, t1 + 5.0);

pt_t2  = pt_fun(t2);
vpt_t2 = vpt_fun(t2);
R1 = norm(pt_t2 - pr_t1);

%% ---------- 求 t3 ----------
fun_t3 = @(t) norm(pr_fun(t) - pt_t2) - c0 * (t - t2);

t3_init = t2 + R1 / c0;
t3 = local_fzero_with_expand(fun_t3, t3_init, t2, t2 + 5.0);

pr_t3  = pr_fun(t3);
vpr_t3 = vpr_fun(t3);
R2 = norm(pr_t3 - pt_t2);

%% 双程斜距
R_two_way = R1 + R2;

%% 非停走多普勒
% 按你给的公式对应写法：
% fd = 1/lambda * [ u1·(vr(t1)-vt(t2)) + u2·(vt(t2)-vr(t3)) ]
u1 = (pt_t2 - pr_t1) / R1;
u2 = (pr_t3 - pt_t2) / R2;

fd_val = ( dot(u1, (vpr_t1 - vpt_t2)) + dot(u2, (vpt_t2 - vpr_t3)) ) / lambda;

end


%% ============================================================
% 带区间扩展的 fzero 求根
%% ============================================================
function root = local_fzero_with_expand(fun, x0, left_min, right_max)

try
    root = fzero(fun, x0);
    if isreal(root) && root >= left_min && root <= right_max
        return;
    end
catch
end

step_list = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0];

for s = step_list
    a = max(left_min, x0 - s);
    b = min(right_max, x0 + s);

    fa = fun(a);
    fb = fun(b);

    if ~isnan(fa) && ~isnan(fb) && isfinite(fa) && isfinite(fb)
        if fa == 0
            root = a;
            return;
        elseif fb == 0
            root = b;
            return;
        elseif fa * fb < 0
            root = fzero(fun, [a, b]);
            return;
        end
    end
end

Ngrid = 2000;
grid = linspace(left_min, right_max, Ngrid);
vals = nan(size(grid));

for i = 1:Ngrid
    vals(i) = fun(grid(i));
end

idx = find(isfinite(vals(1:end-1)) & isfinite(vals(2:end)) & ...
           vals(1:end-1).*vals(2:end) <= 0, 1);

if ~isempty(idx)
    root = fzero(fun, [grid(idx), grid(idx+1)]);
    return;
end

error('fzero 求根失败：在给定区间内没有找到可靠根。');
end
