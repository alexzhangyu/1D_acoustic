function [rl,lhsl] = elm_mtx_solid(yaf_1, yaf_2,...
                                   yam_1, yam_2,...
                                   aa,   bb,...
                                   rho_0, cv_0,... 
                                   k,   sh_mod,...
                                   h,...
                                   b11,...
                                   b11_dot,...
                                   itau_s)
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
global afb
global amb
global gamab
global c_s
%% some fundamental matrix
M = h.*[1.0/3.0 1.0/6.0; 1.0/6.0 1.0/3.0];
E = [-0.5 -0.5; 0.5 0.5];
P = (1.0/h).*[1.0 -1.0; -1.0 1.0];
Q = [-1.0;1.0];
%%
% A0
A0 = zeros(3,3);
A0(1,1) = rho_0*bb;
A0(1,3) = rho_0*aa;
A0(2,2) = rho_0;
A0(3,3) = rho_0*cv_0;
% A1
A1 = zeros(3,3);
A1(1,2) = rho_0;
A1(2,1) = 1.0;
A1(3,2) = p0;
%
% A1^T
A1T = transpose(A1);
%
K = zeros(3,3);
K(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
K(3,3) = k;
%
L = zeros(3,1);
L(2,1) = (2.0/3.0)*sh_mod*( b11 + (dt*afb*(amb-gamab))/amb*b11_dot);
%
% Tau
if(itau_s == -1)
    Tau = zeros(3,3);
elseif(itau_s == 1)    
    Tau(1,1) = h/2.0*c_s;
    Tau(2,2) = h/(2.0*c_s*rho_0);
    Tau(3,3) = h/(2.0*c_s*rho_0*cv_0);
else
    disp('Error: set correct Tau for solid')
end
%
rl = zeros(nshl*nflow,1);
lhsl = zeros(nshl*nflow,nshl*nflow);
%
for irow = 1:nshl
   rl( nflow*(irow-1)+1:nflow*irow,1) = ...
       (M(irow,1).*(A0*yam_1)...
      + M(irow,2).*(A0*yam_2))...
      -( E(irow,1).*(A1*yaf_1)...
        +E(irow,2).*(A1*yaf_2) )...
      +( P(irow,1).*(K*yaf_1)...
        +P(irow,2).*(K*yaf_2) )...
      + Q(irow,1).*L;
end
%adding the stabilization term to res
for irow = 1:nshl
    rl( nflow*(irow-1)+1:nflow*irow,1) = ...
      rl( nflow*(irow-1)+1:nflow*irow,1)...
     +(  E(irow,1)*(A1*Tau*A0*yam_1)...
       + E(irow,2)*(A1*Tau*A0*yam_2) )...
     +(  P(irow,1)*(A1*Tau*A1*yaf_1)...
        +P(irow,2)*(A1*Tau*A1*yaf_2) );  
%      +(  P(irow,1)*(A1*Tau*A1*yaf_1)...
%         +P(irow,2)*(A1*Tau*A1*yaf_2) );
%          +(  E(irow,1)*(A1*Tau*A0*yam_1)...
%        + E(irow,2)*(A1*Tau*A0*yam_2) )...
%      +(  P(irow,1)*(A1*Tau*A1*yaf_1)...
%         +P(irow,2)*(A1*Tau*A1*yaf_2) );  
        
end
%
for irow = 1:nshl
    for jclm = 1:nshl
        itemp = nflow*(irow-1)+1;
        jtemp = nflow*(jclm-1)+1;
        lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
            (am/gama/dt/af)*M(irow,jclm).*A0...
           - E(irow,jclm).*A1...
           + P(irow,jclm).*K;
            
    end
end
%
% adding the stabilization to LHS
for irow = 1:nshl
    for jclm = 1:nshl
        itemp = nflow*(irow-1)+1;
        jtemp = nflow*(jclm-1)+1;
        lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
           lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) )...
           + P(irow,jclm).*(A1*Tau*A1);
%        lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) ) =...
%            lhsl( itemp:(itemp+nflow-1),jtemp:(jtemp+nflow-1) )...
%            + (am/(gama*dt*af))*E(irow,jclm).*(A1*Tau*A0)...
%            + P(irow,jclm).*(A1*Tau*A1);
    end
end
   
end