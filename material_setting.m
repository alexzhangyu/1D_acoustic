function material_setting()
global p0
global T0
% fluid %assuming air
global mu_l
global k_l
%% NOTICE,ONLY USED FOR TEST CASE WATER+AIR+DG
global mu_w
global k_w
mu_w = 1.0e-3; %(Pa*s)
k_w = 0.6; %water (w/(K m)) 
%%
mu_l = 1.8e-5;%1.8e-5;
%   mu_l = 1.8;
% mu_l = 1.0e-3; %water
k_l = 0.026;
% k_l = 0.6; %water
R = 8.314e3;
mv = 29; %molecular weight for air g/mol
rr = 1.4;
% solid
global bk_mod
global sh_mod
global E_mod
global k_s
bk_mod = 10.0e6;
sh_mod = 15.0e6;
E_mod = 2.0*sh_mod;
k_s = 31.1;
% common propterties: aa:related to temp
%                     bb: related to pressure
%                     rho_0 : intial density
%                     cv : specific heat per vol
% in the common properties arrays, 1st entry: fluid
%                                  2nd entry: solid
%                                  3rd entry: water 
global rho_0
global aa
global bb
global cv_0
rho_0 = zeros(3,1);
aa = zeros(3,1);
bb = zeros(3,1);
cv_0 = zeros(3,1);
% 
rho_0(1,1) = 1.2111;
% rho_0(1,1) = 1000; %water
rho_0(2,1) = 1.8775e3;
%%NOTICE, ONLY FOR CASE WITH WATER+AIR+DG
rho_0(3,1) = 100;%water
% 
aa(1,1) = -1.0/T0;
% aa(1,1) = -1.0e-4;
aa(2,1) = -3.0*2.22e-5;
%%NOTICE, ONLY FOR CASE WITH WATER+AIR+DG
aa(3,1) = -1.0e-4;
% 
bb(1,1) = 1.0/p0;
% bb(1,1) = 4e-10;
bb(2,1) = 1.0/bk_mod;
%%NOTICE, ONLY FOR CASE WITH WATER+AIR+DG
bb(3,1) = 1e-7;
%
cv_0(1,1) = R/mv/(rr-1.0);
% cv_0(1,1) = 4100;
cv_0(2,1) = 1586;
%%NOTICE, ONLY FOR CASE WITH WATER+AIR+DG
cv_0(3,1) = 1460;

% longitude wave speed for solid
global c_s
 c_s = sqrt(E_mod/rho_0(2,1));
% c_s = sqrt(bk_mod/rho_0(2,1));
%%NOTICE, ONLY FOR CASE WITH WATER+AIR+DG
global c_w
 c_w = sqrt(1.0/bb(3,1)/rho_0(3,1));
% 
end