function [ K ] = setK_fluid(K, mu,lamda )
%fill the element K matrix
%   Using the proper material properties
%
K(2,2) = (4.0/3.0)*mu;
K(3,3) = lamda;
% % %% debugging
% K(3,:) = 0.0;
% K(:,3) = 0.0;
end