function [rl,lhsl] = elm_mtx_midpoint( y1, y2,...
                                       aa,    bb,...
                                       rho_0, cv_0,... 
                                       k,     mu,...
                                       h)
%%
%calculating the element level residual and lhs matrix
%%
global nflow
global nshl
global p0
global af
global am
global gama
global dt
global c_air
%% some fundamental matrix
M = h*[1.0/3.0 1.0/6.0; 1.0/6.0 1.0/3.0];
E = [-0.5 -0.5; 0.5 0.5];
P = (1.0/h)*[1.0 -1.0; -1.0 1.0];
%%
% A0
% A0 = zeros(3,3);
% A0(1,1) = rho_0*bb;
% A0(1,3) = rho_0*aa;
% A0(2,2) = rho_0;
% A0(3,3) = rho_0*cv_0;
A0 = eye(3);
% A1
A1 = zeros(3,3);
% A1(1,1) = 0.5;
A1(1,2)= 0.5;
A1(2,1)= 0.5;
% A1(1,2) = rho_0;
% A1(2,1) = 1.0;
% A1(3,2) = p0;
% A1^T
A1T = transpose(A1);
%
K = zeros(3,3);
K(1,1) = 1e-4;
K(2,2) = 1e-4;
% K(2,2) = (4.0/3.0)*mu;
% K(3,3) = k;
% K(1,1) = A1(1,1)*h/8.0;
% Tau
Tau = zeros(3,3);
% Tau(1,1) = h/(2.0*0.5);
% Tau(2,2) = h/(2.0*0.5);
% Tau(3,3) = h/(2.0*0.5);
% Tau(1,1) = h/(2.0*c_air);
% Tau(2,2) = h/(2.0*c_air);
% Tau(3,3) = h/(2.0*c_air);
%
rl = zeros(nshl*nflow,1);
lhsl = zeros(nshl*nflow,nshl*nflow);
%
for irow = 1:nshl
   rl( nflow*(irow-1)+1:nflow*irow,1) = ...
        2.0*M(irow,1)*(A0*y1)...
      + 2.0*M(irow,2)*(A0*y2)...
      + dt*( E(irow,1)*(A1*y1)...
        +E(irow,2)*(A1*y2) )...
      - dt*( P(irow,1)*(K*y1)...
        +P(irow,2)*(K*y2) );
end
% bk euler
% for irow = 1:nshl
%    rl( nflow*(irow-1)+1:nflow*irow,1) = ...
%         M(irow,1)*(A0*y1)...
%       + M(irow,2)*(A0*y2);
% end
%
% applying the stabilization
for irow = 1:nshl
   rl( nflow*(irow-1)+1:nflow*irow,1) = ...
      rl( nflow*(irow-1)+1:nflow*irow,1)...
      + 2.0*E(irow,1)*(A1T*Tau*A0*y1)...
      + 2.0*E(irow,2)*(A1T*Tau*A0*y2)...
      - dt*( P(irow,1)*(A1T*Tau*A1*y1)...
            +P(irow,2)*(A1T*Tau*A1*y2) );
end
%%LHS
for irow = 1:nshl
    for jclm = 1:nshl
        itemp = nflow*(irow-1)+1;
        jtemp = nflow*(jclm-1)+1;
        lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
            2.0*M(irow,jclm)*A0...
           - dt*E(irow,jclm)*A1...
           + dt*P(irow,jclm)*K;
            
    end
end
% bk euler
% for irow = 1:nshl
%     for jclm = 1:nshl
%         itemp = nflow*(irow-1)+1;
%         jtemp = nflow*(jclm-1)+1;
%         lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
%              M(irow,jclm)*A0...
%            - dt*E(irow,jclm)*A1...
%            + dt*P(irow,jclm)*K;
%             
%     end
% end
% applying the stabilization
for irow = 1:nshl
    for jclm = 1:nshl
        itemp = nflow*(irow-1)+1;
        jtemp = nflow*(jclm-1)+1;
        lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
          lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) )...
           + 2.0*E(irow,jclm)*(A1T*Tau*A0)...
           + dt*P(irow,jclm)*(A1T*Tau*A1);
            
    end
end

   
end