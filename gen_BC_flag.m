function [ibcg,ibcb,bcg,bcb] = gen_BC_flag()
% In the 1D case, the boundary conditions only needed to assigned at two
% end nodes or dof(1:nflow,1) and dof(Nnode*nflow-2:Nnode*nflow,1) in terms
% of global dof
% ibcg(1,1:3) : dirichlet BC flag of p,v,T of left end, 0-not set,1-set
% ibcg(2,1:3) : dirichlet BC on p,v,T of right end, 0-not set, 1-set
% 
% bcg(1,1:3): dirichlet BC of p,v,T of left end
% bcg(2,1:3): dirichlet BC of p,v,T of right end
% 
% ibcb(1,1:3) : Neumann BC flag on p,v,T of left end, 0- not set, 1-set
% ibcb(1,1:3) : Neumann BC flag on p,v,T of right end, 0- not set, 1-set
% 
% bcb(1,1:3): Neumann BC of p,v,T of left end
% bcb(2,1:3): Neumann BC of p,v,T of right end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
global nflow
global c_air
global omg
global dt
global bk_mod
global E_mod
ibcg = uint8(zeros(2,nflow));
bcg = zeros(2,nflow);
ibcb = uint8(zeros(2,nflow));
bcb = zeros(2,nflow);
%% dirichlet BC  
%% left end
   ibcg(1,1) = 1;
   ibcg(1,2) = 1;
   ibcg(1,3) = 1;
% BC values
%    bcg(1,1) = 4.117714504685213e+02;
   bcg(1,1) = cal_eigenvalue;
%    bcg(1,1) = 0.002873492036850;
%    bcg(1,1)= 0.0;
%    bcg(1,1) = 1.370220492929288e+05; %sqrt(bk_mod/E_mod)
   bcg(1,2)= 0.1;
   bcg(1,3) = 0.0;
%    bcg(1,3)=0.338837709155365;
%% right end   
   ibcg(2,1) = 0;
   ibcg(2,2) = 1;
   ibcg(2,3) = 1;
% BC values
   bcg(2,1) = 0.0;
   bcg(2,2) = 0.0;
   bcg(2,3) = 0.0;
%    bcg(2,3) = bcg(1,3);
%% Neumann BC
%% left end
   ibcb(1,1) = 0;
   ibcb(1,2) = 0;
   ibcb(1,3) = 0;
% % BCB values
%    bcb(1,1) = 1.0;
%    bcb(1,2) = 1.0;
%    bcb(1,3) = 0.0;
%% right end   
   ibcb(2,1) = 0;
   ibcb(2,2) = 0;
   ibcb(2,3) = 0;
% % BC values
%    bcb(2,1) = 0.0;
%    bcb(2,2) = 0.0;
%    bcb(2,3) = 0.0;

  
end