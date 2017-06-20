function [A0,A1,K] = set_elmmtx_water(A0,A1,K)
% setting the A1,K and L for the element which is water
%
global p0
global mu_w
global k_w
global rho_0
global aa
global bb
global cv_0
%
% A0
% A0(1,1) = rho_0(3,1)*bb(3,1);
% A0(1,3) = rho_0(3,1)*aa(3,1);
% A0(2,2) = rho_0(3,1);
% A0(3,3) = rho_0(3,1)*cv_0(3,1);
A0 = setA0(A0, rho_0(3,1),aa(3,1),bb(3,1),cv_0(3,1) );
%  A1
% A1(1,2) = rho_0(3,1);
% A1(2,1) = 1.0;
% A1(3,2) = p0;
A1 = setA1(A1, rho_0(3,1), p0);
% %%debugging
% A1(3,:) = 0.0;
% A1(:,3) = 0.0;
% K
% K(2,2) = (4.0/3.0)*mu_w;
% K(3,3) = k_w;
K = setK_fluid(K, mu_w,k_w );
% %%debugging
% K(3,:) = 0.0;
% K(:,3) = 0.0;
end