function [Y] = impose_BC_midpoint(Y,itime)
%% enforcing the BC
global nflow
global Nnode
global dt
global omg
%
global ibcg
global bcg
% left end
for i = 1:nflow
    if(ibcg(1,i)==1)
        Y(i,1) = bcg(1,i);
    end
end
% % pressure on left
% %     Y(1,1) = 1.0*sin(omg*itime*dt);
%     Y(1,1) = 1.0;
% %     Y(2,1) = 1.0*sin(omg*itime*dt);
%     Y(2,1) = 1.0;
%     Y(3,1) = 0.0;
% fixed on the right
for i = 1:nflow
    if(ibcg(2,i)== 1)
        Y(nflow*Nnode-3+i,1) = bcg(2,i);
    end
end
%     Y(nflow*Nnode,1) = 0.0;% T on the right
%      Y(nflow*Nnode-1,1) = 0.0;%v on the right
%     Y(nflow*Nnode-2,1) = 0.0; % p on the right
end