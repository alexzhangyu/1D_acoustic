close all
clc
clear
global p0
global T0
p0=1e5;
T0=288;
%
material_setting
% fluid %assuming air
global mu_l
global k_l
% solid
global bk_mod
global sh_mod
global E_mod
global k_s
global c_s
%
global rho_0
global aa
global bb
global cv_0
%
afb = 1.0;
gamab = 1.0;
amb = 1.0;
%
CFL = 2.0;
h = 1e-2;
dt = h/c_s*CFL;
%
material_setting
%
% A0
A0 = zeros(3,3);
A0(1,1) = rho_0(2,1)*bb(2,1);
A0(1,3) = rho_0(2,1)*aa(2,1);
A0(2,2) = rho_0(2,1);
A0(3,3) = rho_0(2,1)*cv_0(2,1);
% A0(:,3) = 0.0;
% A0(3,:) = 0.0;
% A1
A1 = zeros(3,3);
A1(1,2) = rho_0(2,1);
A1(2,1) = 1.0;
A1(3,2) = p0;
%
K = zeros(3,3);
K(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
K(3,3) = k_s;
%
a0a1 = (A0)^-1*A1;
[V,D] = eig(a0a1);
K_z = V^-1*(A0)^-1*K*V;
% vv1 = [0;0;1];
% vv2 = [-(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1))/bb(1,1)/p0; sqrt(-bb(1,1)*cv_0(1,1)*(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1)))/bb(1,1)/p0;1];
% vv3 = [-(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1))/bb(1,1)/p0; -sqrt(-bb(1,1)*cv_0(1,1)*(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1)))/bb(1,1)/p0;1];
%%
syms y1 y2 y3
z=V^-1*[y1;y2;y3];
f=subs(z(3),[y2,y3],[1,0]);
y_fill = solve(f,y1);
disp(eval(y_fill));
%%
% syms y1 y2
% z=V'*[y1;y2];
%  f=subs(z(2),y2,1);
%  y_fill = solve(f,y1);
%  disp(eval(y_fill));