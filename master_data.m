clc
clear
close all

% load opt_all.mat

d_star          = 3*.249e-9; 
cell_res        = 5e-6;
cell_max        = (sqrt(2)/2)*cell_res*2;
cell_maxd       = cell_max/2;
disl_num        = 1;
% bvec            = [1,0,0];
% t               = 0;

res_r           = 1000000;
res_a           = 1000000;
res_t           = 1000000;
res_d           = 1000000;

master_sigma; 
fun_norm        = @(x)             (sum(x.^2,2)).^0.5;
fun_scale       = @(x,lb,ub)       ((x-(ones(size(x,1),1)*min(x)))./((ones(size(x,1),1)*max(x))-(ones(size(x,1),1)*min(x))))*(ub-lb)+lb; 
fun_scaleback   = @(x,a,b,I,M)     ((x-a)./(b-a)).*(ones(size(x,1),1)*M-ones(size(x,1),1)*I)+ones(size(x,1),1)*I;
fun_sigma       = @(sxx,sxy,syy,b) b(:,1).*(sxx.*conj(b(:,1)) + sxy.*conj(b(:,2))) + b(:,2).*(sxy.*conj(b(:,1)) + syy.*conj(b(:,2)));
fun_fg          = @(a,r,bx,by,t)   -bx.*(conj(by).*(((cos(t).^2 - sin(t).^2).*(cos(a) - 2.*cos(a).^3))./r - cos(t).*sin(t).*((3.*sin(a) - 2.*sin(a).^3)./r + (sin(a) - 2.*sin(a).^3)./r)) + conj(bx).*((cos(t).^2.*(3.*sin(a) - 2.*sin(a).^3))./r - (sin(t).^2.*(sin(a) - 2.*sin(a).^3))./r + (2.*cos(t).*sin(t).*(cos(a) - 2.*cos(a).^3))./r)) - by.*(conj(bx).*(((cos(t).^2 - sin(t).^2).*(cos(a) - 2.*cos(a).^3))./r - cos(t).*sin(t).*((3.*sin(a) - 2.*sin(a).^3)./r + (sin(a) - 2.*sin(a).^3)./r)) + conj(by).*((sin(t).^2.*(3.*sin(a) - 2.*sin(a).^3))./r - (cos(t).^2.*(sin(a) - 2.*sin(a).^3))./r + (2.*cos(t).*sin(t).*(cos(a) - 2.*cos(a).^3))./r));
fun_dist        = @(a1,a2,r1,r2)   sqrt( r1.^2 + r2.^2 - 2.*r1.*r2.*cos(a2-a1) );

%% compute the inverse of the stress functions 
syms r a t
sxx_inv         = finverse(sigma.polar.global.xx(r,a,t),r);
sxy_inv         = finverse(sigma.polar.global.xy(r,a,t),r);
syy_inv         = finverse(sigma.polar.global.yy(r,a,t),r);

%% compute the max and min of the function 
% extereme points at d_star
e1 = simplify(diff(sigma.polar.global.xx(d_star,a,t),a));
e2 = simplify(diff(sigma.polar.global.xx(d_star,a,t),t));
q  = solve(e1==0,e2==0);
p1 = [q.a,q.t];

e1 = simplify(diff(sigma.polar.global.xy(d_star,a,t),a));
e2 = simplify(diff(sigma.polar.global.xy(d_star,a,t),t));
q  = solve(e1==0,e2==0);
p2 = [q.a,q.t];

e1 = simplify(diff(sigma.polar.global.yy(d_star,a,t),a));
e2 = simplify(diff(sigma.polar.global.yy(d_star,a,t),t));
q  = solve(e1==0,e2==0);
p3 = [q.a,q.t];

% all extreme points at d_start (we are interested in maximum absolute of
% sigmas for these points 
p_d_star  = unique([eval(p1);eval(p2);eval(p3)],'rows');

% extreme points at cell_max (radius of the cell)
e1 = simplify(diff(sigma.polar.global.xx(cell_max,a,t),a));
e2 = simplify(diff(sigma.polar.global.xx(cell_max,a,t),t));
q  = solve(e1==0,e2==0);
p1 = [q.a,q.t];

e1 = simplify(diff(sigma.polar.global.xy(cell_max,a,t),a));
e2 = simplify(diff(sigma.polar.global.xy(cell_max,a,t),t));
q  = solve(e1==0,e2==0);
p2 = [q.a,q.t];

e1 = simplify(diff(sigma.polar.global.yy(cell_max,a,t),a));
e2 = simplify(diff(sigma.polar.global.yy(cell_max,a,t),t));
q  = solve(e1==0,e2==0);
p3 = [q.a,q.t];

% all extreme points at cell_max (we are interested in minimum absolute of
% sigmas for these points 
p_cell_max  = unique([eval(p1);eval(p2);eval(p3)],'rows');

% compute the minimum and maximum values of sigmas within the simulation
% cell by investigating all the extreme points 

s_lim_xx        = [min(abs(sigma.polar.global.xx(cell_max,p_cell_max(:,1),p_cell_max(:,2)))),max(abs(sigma.polar.global.xx(d_star,p_d_star(:,1),p_d_star(:,2))))];
s_lim_xy        = [min(abs(sigma.polar.global.xy(cell_max,p_cell_max(:,1),p_cell_max(:,2)))),max(abs(sigma.polar.global.xy(d_star,p_d_star(:,1),p_d_star(:,2))))];
s_lim_yy        = [min(abs(sigma.polar.global.yy(cell_max,p_cell_max(:,1),p_cell_max(:,2)))),max(abs(sigma.polar.global.yy(d_star,p_d_star(:,1),p_d_star(:,2))))];

%% generate r and a randomely such the functions be uniform and between the minimum and maximum values of that function (between d_star and cell_max); for 
% xx
s               = rng;
rng(s);
t_xx            = rand(res_t,disl_num)*2*pi;
a_xx            = rand(res_a,disl_num)*2*pi;
r_xx            = ((rand(res_r,disl_num) > 0.5)*2 - 1).*fun_scale(rand(res_r,disl_num),s_lim_xx(1),s_lim_xx(2));

% xy
s               = rng;
rng(s);
t_xy            = rand(res_t,disl_num)*2*pi;
a_xy            = rand(res_a,disl_num)*2*pi;
r_xy            = ((rand(res_r,disl_num) > 0.5)*2 - 1).*fun_scale(rand(res_r,disl_num),s_lim_xy(1),s_lim_xy(2));

% yy
s               = rng;
rng(s);
t_yy            = rand(res_t,disl_num)*2*pi;
a_yy            = rand(res_a,disl_num)*2*pi;
r_yy            = ((rand(res_r,disl_num) > 0.5)*2 - 1).*fun_scale(rand(res_r,disl_num),s_lim_yy(1),s_lim_yy(2));

%% Generate the radiuses 
r = r_xx; a = a_xx; t = t_xx;
r_xx            = eval(sxx_inv);%sigma.polar.global.xx(r_xx,a_xx,t_xx);

r = r_xy; a = a_xy; t = t_xy;
r_xy            = eval(sxy_inv);%sigma.polar.global.xy(r_xy,a_xy,t_xy);

r = r_yy; a = a_yy; t = t_yy;
r_yy            = eval(syy_inv);%sigma.polar.global.yy(r_yy,a_yy,t_yy);

%% phase modification
% turn every negetive r to positive and change the corresponding a
ind_neg = find(r_xx<0); r_xx(ind_neg) = -r_xx(ind_neg); a_xx(ind_neg) = wrapTo2Pi(a_xx(ind_neg)+pi);
ind_neg = find(r_xy<0); r_xy(ind_neg) = -r_xy(ind_neg); a_xy(ind_neg) = wrapTo2Pi(a_xy(ind_neg)+pi);
ind_neg = find(r_yy<0); r_yy(ind_neg) = -r_yy(ind_neg); a_yy(ind_neg) = wrapTo2Pi(a_yy(ind_neg)+pi);

%% remove radiuses less than d_star 
ind_xx = find(abs(r_xx)<d_star | abs(r_xx)>cell_max);
ind_xy = find(abs(r_xy)<d_star | abs(r_xy)>cell_max);
ind_yy = find(abs(r_yy)<d_star | abs(r_yy)>cell_max);

ind_bad = unique([ind_xx;ind_xy;ind_yy]);

r_xx(ind_bad,:) = [];
r_xy(ind_bad,:) = [];
r_yy(ind_bad,:) = [];

a_xx(ind_bad,:) = [];
a_xy(ind_bad,:) = [];
a_yy(ind_bad,:) = [];

t_xx(ind_bad,:) = [];
t_xy(ind_bad,:) = [];
t_yy(ind_bad,:) = [];

r      = [r_xx;r_xy;r_yy];
a      = [a_xx;a_xy;a_yy];
t      = [t_xx;t_xy;t_yy];

sxx    = sigma.polar.global.xx(r,a,t);
sxy    = sigma.polar.global.xy(r,a,t);
syy    = sigma.polar.global.yy(r,a,t);

s      = [sxx,sxy,syy];

%master_plots

%% Inputs & Outputs 
r_log  = log10(r);
r_max  = log10(d_star);
a_max  = 2*pi;
t_max  = 2*pi;
smaxoff=[];

s_max  = max(abs(s(:)));

datain_raw       = [a,r_log,t];
dataou_raw       = [s];

datain_scl       = [a./a_max,r_log./r_max,t./t_max];
dataou_scl       = [s./s_max];

%% save 
save('data_raw','datain_raw','dataou_raw','r_max','a_max','t_max','s_max','-v7.3')
save('data_scl','datain_scl','dataou_scl','r_max','a_max','t_max','s_max','-v7.3')

% save in dict in .txt for python
% raw
datain_list = fun_mat2list(datain_raw);
dataou_list = fun_mat2list(dataou_raw);
dict_val    = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid         = fopen('data_inou_raw.txt','wt');
fprintf(fid, dict_val);
fclose(fid);

% scl
datain_list = fun_mat2list(datain_scl);
dataou_list = fun_mat2list(dataou_scl);
dict_val    = append("{'datain': ",string(datain_list),", 'dataou': ", string(dataou_list), "}");
fid         = fopen('data_inou_scl.txt','wt');
fprintf(fid, dict_val);
fclose(fid);



