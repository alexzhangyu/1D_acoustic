function [y_fill] = cal_eigenvalue
clc
clear
%
global p0
global T0
p0=1e5;
T0=288;
% fluid %assuming air
global mu_l
global k_l
% solid
global bk_mod
global sh_mod
global E_mod
global k_s
%
global rho_0
global aa
global bb
global cv_0
%
material_setting
%
% A0
A0 = zeros(3,3);
A0(1,1) = rho_0(1,1)*bb(1,1);
A0(1,3) = rho_0(1,1)*aa(1,1);
A0(2,2) = rho_0(1,1);
A0(3,3) = rho_0(1,1)*cv_0(1,1);
% debugging, pure advective
% A0(1,3) = 0.0;
% A0(3,3) = 1.0;
% A1
A1 = zeros(3,3);
A1(1,2) = rho_0(1,1);
A1(2,1) = 1.0;
A1(3,2) = p0;
% debugging, pure advective
%  A1(:,3) = 0.0;
%  A1(3,:) = 0.0;
%
a0a1 = (A0)^-1*A1;
[V,D] = eig(a0a1);
% vv1 = [0;0;1];
% vv2 = [-(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1))/bb(1,1)/p0; sqrt(-bb(1,1)*cv_0(1,1)*(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1)))/bb(1,1)/p0;1];
% vv3 = [-(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1))/bb(1,1)/p0; -sqrt(-bb(1,1)*cv_0(1,1)*(aa(1,1)*p0-cv_0(1,1)*rho_0(1,1)))/bb(1,1)/p0;1];
%%
syms y1 y2 y3
z=V^-1*[y1;y2;y3];
f=subs(z(3),[y2,y3],[0.1,0]);
y_fill = solve(f,y1);
disp(eval(y_fill));
%%
% syms y1 y2
% z=V'*[y1;y2];
%  f=subs(z(2),y2,1);
%  y_fill = solve(f,y1);
%  disp(eval(y_fill));
end
