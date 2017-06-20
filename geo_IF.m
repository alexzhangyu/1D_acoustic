function [IEN,IENIF0,IENIF1,x,h,Nel,Nnode] = geo_IF(x0, x1,xif, N0,N1)
%
  h0 = ( xif - x0)/N0; % mesh size in phase0
  h1 = ( x1 - xif)/N1; % mesh size in phase1
  h = zeros(2,1); % h(1,1) for fluid
                 % h(2,1) for solid
  h(1,1) = h1; % fluid on the right
  h(2,1) = h0; % solid on the left
%
  Nel = N0 + N1;
  Nnode = N0 + N1 + 2; %accouting the duplicated node 
                     % at interface
%
% global IEN array
  IEN = zeros(Nel,3);
  for i = 1:N0
    IEN(i,1) = i;
    IEN(i,2) = i+1;
    IEN(i,3) = 0; %solid phase on the left
  end
  for i = (N0+1):Nel
    IEN(i,1)= i+1;
    IEN(i,2) = i+2;
    IEN(i,3) = 1; %fluid phase
  end
% IENIF array
    IENIF0 = zeros(1,3);
    IENIF1 = zeros(1,3);
%
    IENIF0(1,1) = N0; % 0 side, close to the interface
    IENIF0(1,2) = N0+1; % 0 side, on the interface
    IENIF0(1,3) = 0; % 0 side is solid
    IENIF1(1,1) = N0+2; % 1 side, on the interface
    IENIF1(1,2) = N0+3; % 1 side, close to the interface
    IENIF1(1,3) = 1; % 1 side is fluid
% global coords
  x = zeros(Nnode,1);
% filling global coords
  for i = 1:N0+1
    x(i) = x0+(i-1)*h0;
  end
  for i= N0+2:Nnode
    x(i) = xif + (i-(N0+2))*h1;
  end    
end