function [ L ] = setL_solid(L, sh_mod,B11,B11_dot )
%fill the element L matrix
%   Using the proper material properties
%%
global dt
global afb
global amb
global gamab
%
L(2,1) = (2.0/3.0)*sh_mod*( B11 + (dt*afb*(amb-gamab))/amb*B11_dot);
end