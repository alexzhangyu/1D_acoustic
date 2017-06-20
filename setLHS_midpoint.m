function [LHS] = setLHS_midpoint(LHS)
%%%%%%%%%%%%%
global nflow
global Nnode
%
global ibcg
global bcg
% setting p,u,t on the left
for i = 1:nflow
    if(ibcg(1,i)== 1)
        LHS(i,:) = 0.0;
        LHS(i,i) = 1.0;
    end
end
%    LHS(1,:) = 0.0;
%    LHS(2,:) = 0.0;
%    LHS(3,:) = 0.0;
% %   
%      LHS(1,1) = 1.0;
% %     
%      LHS(2,2) = 1.0;
% %     
%      LHS(3,3) = 1.0;
% setting the u and t = 0 on the right 
for i = 1:nflow
    if(ibcg(2,i)== 1)
        LHS((nflow*Nnode-3+i),:) = 0.0;
        LHS((nflow*Nnode-3+i),(nflow*Nnode-3+i)) = 1.0;
    end
end
%    LHS(nflow*Nnode,:) = 0.0; %T
%    LHS((nflow*Nnode-1),:) = 0.0; %v
%    LHS((nflow*Nnode-2),:) = 0.0; %p
% % % T RHS
%    LHS(nflow*Nnode,nflow*Nnode) = 1.0; 
% % % v RHS  
%    LHS(nflow*Nnode-1,nflow*Nnode-1) = 1.0;
% % % p RHS   
%    LHS(nflow*Nnode-2,nflow*Nnode-2) = 1.0;

end