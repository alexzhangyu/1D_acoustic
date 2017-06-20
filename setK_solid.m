function [ K ] = setK_solid(K, sh_mod,lamda )
%fill the element K matrix
%   Using the proper material properties
%%
global dt
global afb
global amb
global gamab
%
K(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
K(3,3) = lamda;
end