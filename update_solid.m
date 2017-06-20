function [b11,b11_dot] = update_solid( y_af,...
                            N_solid,ien,...
                            b11,  b11_dot,...
                            h)
% Update the left cauchy green tensor b for solid elements
%%
global nflow
%
global dt
global afb
global amb
global gamab
%
    for iel = 1:N_solid
       up1s = (ien(iel,1)-1)*nflow+1; % starting location for node1
       up1e = ien(iel,1)*nflow; % end location for node1
       up2s = (ien(iel,2)-1)*nflow+1; % starting location for node2
       up2e = ien(iel,2)*nflow; %end location for node2
% localized Y       
       Yaf_1s = y_af( up1s : up1e,1);
       Yaf_2s = y_af( up2s : up2e,1 );
% linearized approxi for grad Y
       gradv_af = (Yaf_2s(2,1)-Yaf_1s(2,1))/h;
% update b11 and b11_dot
       b11(iel,1) = b11(iel,1) + dt*(amb-gamab)/amb*b11_dot(iel,1)...
                 + 2.0*(gamab*dt/amb)*gradv_af;
       b11_dot(iel,1) = (amb-1.0)/amb*b11_dot(iel,1)...
                    + (2.0/amb)*gradv_af;
 
    end


end

