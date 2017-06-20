function [ A1 ] = setA1(A1, rho_0,p )
%fill the element A1 matrix
%   Using the proper material properties
%
A1(1,2) = rho_0;
A1(2,1) = 1.0;
A1(3,2) = p;
% %% debugging
% A1(1,2) = 1;
% A1(2,1) = 1;
% A1(3,:) = 0.0;
% A1(:,3) = 0.0;
end

