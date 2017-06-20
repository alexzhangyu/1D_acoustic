function [IEN,x,h,Nel,Nnode] = geo_CG(x0, x1,N)
%
  h = zeros(2,1); % h(1,1) for fluid
                 % h(2,1) for solid
  h(1,1) = ( x1 - x0 )/N; % fluid, default
  h(2,1) = ( x1 - x0 )/N; % solid, optional
%
  Nel = N;
  Nnode = N+1;     % total number of nodes
%
% global IEN array
  IEN = zeros(Nel,3);
  for i = 1:N
    IEN(i,1) = i;
    IEN(i,2) = i+1;
    IEN(i,3) = 0; % fluid, default !changed to solid
  end
% global coords
  x = zeros(Nnode,1);
% filling global coords
  for i = 1:Nnode
    x(i,1) = x0+(i-1)*h(1,1);
  end
  
end