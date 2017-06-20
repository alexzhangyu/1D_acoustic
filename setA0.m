function [ A0 ] = setA0( A0, rho_0,aa,bb,cv_0  )
%fill the element A0 matrix
%   Using the proper material properties
%%
% A0
A0(1,1) = rho_0*bb;
A0(1,3) = rho_0*aa;
A0(2,2) = rho_0;
A0(3,3) = rho_0*cv_0;
% debugging
% A0 = eye(3,3);
% A0(1,3) = 0.0;
% A0(3,:) = 0.0;
end

