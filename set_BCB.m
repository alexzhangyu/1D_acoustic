function [resl,resr,LHSl,LHSr] = set_BCB(y_l1,y_l2,y_r1,y_r2,...
                                         mater_l,mater_r,...
                                         b11_0, b11_dot_0,...
                                         b11_1, b11_dot_1,...
                                         h_l,   h_r)
% modifying the LHS and RHS if no BC is assigned at the boundary nodes.

global ibcg
global bcg
global ibcb
global bcb
global nflow
global nshl
%
global p0
% global af
% global am
% global gama
global dt
global afb
global amb
global gamab
%
global mu_l
global k_l
global bk_mod
global sh_mod
global k_s
%
global rho_0
%
nl = -1.0; % out normal at left end
nr = 1.0; % out normal at right end 
% setting some elemental matrix
% for the left end
Ml = [1.0 0]; % NbNa at x=left end
Nl = [-0.5 0.5]*(2.0/h_l); %NbNa,1 at x=left end
Ql = 1.0; % Nb at x=left end
% for the right end
Mr = [0 1.0]; %NbNa at x=right end
Nr = [-0.5 0.5]*(2.0/h_r); %NbNa,1 at x=right end
Qr = 1.0; %Nb at x=right end
% forming the A1 and K matrix wrt material type
if(mater_l == 1)
% %A1
%   A1_l = zeros(3,3);
%   A1_l(1,2) = rho_0(2,1);
%   A1_l(2,1) = 1.0;
%   A1_l(3,2) = p0;
% %
%   K_l = zeros(3,3);
%   K_l(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
%   K_l(3,3) = k_s;
% %
%   L_l = zeros(3,1);
%   L_l(2,1) = (2.0/3.0)*sh_mod*( b11_0 + (dt*afb*(amb-gamab))/amb*b11_dot_0);
%% ONLY USED FOR THE TEST CASE WITH WATER+AIR+DG
  A1_l = zeros(3,3);
  K_l = zeros(3,3);
  L_l = zeros(3,1);
  A0_l =zeros(3,1);
%  
  [A0_l,A1_l,K_l] = set_elmmtx_water(A0_l,A1_l,K_l);

elseif(mater_l == 0) 
% A1
  A1_l = zeros(3,3);
%   A1_l(1,2) = rho_0(1,1);
%   A1_l(2,1) = 1.0;
%   A1_l(3,2) = p0;
%     A1_l(1,2) = 1.0;
%     A1_l(2,1) = 1.0;
  A1_l = setA1(A1_l, rho_0(1,1), p0);
%
  K_l = zeros(3,3);
%   K_l(2,2) = (4.0/3.0)*mu_l;
%   K_l(3,3) = k_l;
  K_l = setK_fluid(K_l, mu_l,k_l );
%
  L_l = zeros(3,1);
else
  disp('wrong material for left end')
end
%
if(mater_r == 1)  
% % A1
%   A1_r = zeros(3,3);
%   A1_r(1,2) = rho_0(2,1);
%   A1_r(2,1) = 1.0;
%   A1_r(3,2) = p0;
% %
%   K_r = zeros(3,3);
%   K_r(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
%   K_r(3,3) = k_s;
% %
%   L_r = zeros(3,1);
%   L_r(2,1) = (2.0/3.0)*sh_mod*( b11_1 + (dt*afb*(amb-gamab))/amb*b11_dot_1);
%% ONLY USED FOR THE TEST CASE WITH WATER+AIR+DG
  A1_r = zeros(3,3);
  K_r = zeros(3,3);
  L_r = zeros(3,1);
  A0_r = zeros(3,3);
%  
  [A0_r,A1_r,K_r] = set_elmmtx_water(A0_r,A1_r,K_r);

elseif(mater_r == 0)
% A1
  A1_r = zeros(3,3);
%   A1_r(1,2) = rho_0(1,1);
%   A1_r(2,1) = 1.0;
%   A1_r(3,2) = p0;
%     A1_r(1,2) = 1.0;
%     A1_r(2,1) = 1.0;
  A1_r = setA1(A1_r, rho_0(1,1), p0);
%
  K_r = zeros(3,3);
%   K_r(2,2) = (4.0/3.0)*mu_l;
%   K_r(3,3) = k_l;
  K_r = setK_fluid(K_r, mu_l,k_l );
%   K_r(:,3) = 0.0;
%   K_r(3,:) = 0.0;
%
  L_r = zeros(3,1);
else
  disp('wrong material for right end');
end

%% output array
resl = zeros(nflow,1);
resr = zeros(nflow,1);
%
LHSl = zeros(nflow, nshl*nflow);
LHSr = zeros(nflow, nshl*nflow);
%% local varibles
y_l1_m = zeros(nflow,1);
y_r2_m = zeros(nflow,1);
%% deal with No BC assignment condition
for i=1:nflow
    if((ibcg(1,i)== 1)&&(ibcb(1,i)==1))
        disp('error:cannot assign dirichlet and Neumann at same point');
    elseif((ibcg(1,i)~= 1)&&(ibcb(1,i)~=1))
        y_l1_m(i,1) = y_l1(i,1);
    elseif(ibcb(1,i)==1)
        disp('error: BC type on left not supported');
    end
%
    if((ibcg(2,i)== 1)&&(ibcb(2,i)==1))
        disp('error:cannot assign dirichlet and Neumann at same point');
    elseif((ibcg(2,i)~= 1)&&(ibcb(2,i)~=1))
        y_r2_m(i,1) = y_r2(i,1);
    elseif(ibcb(2,i)==1)
        disp('error: BC type on right not supported');
    end
end
%% local res on left end, already account for solid
resl = (Ml(1,1).*(A1_l*y_l1_m)*nl + Ml(1,2).*(A1_l*y_l2)*nl)...
     - ( Nl(1,1).*(K_l*y_l1_m)*nl + Nl(1,2).*(K_l*y_l2)*nl )...
      - Ql* L_l*nl;
%% local res on right end, already account for solid
resr = (Mr(1,1).*(A1_r*y_r1)*nr + Mr(1,2).*(A1_r*y_r2_m)*nr)...
     - ( Nr(1,1).*(K_r*y_r1)*nr + Nr(1,2).*(K_r*y_r2_m)*nr )...
      - Qr* L_r*nr;
%% local LHS on left end
LHSl(1:nflow,1:nflow) = Ml(1,1)*nl.*A1_l - Nl(1,1)*nl*K_l;
LHSl(1:nflow,nflow+1:nshl*nflow) = Ml(1,2)*nl.*A1_l - Nl(1,2)*nl*K_l;
%% local LHS on right end
LHSr(1:nflow,1:nflow) = Mr(1,1)*nr.*A1_r - Nr(1,1)*nr*K_r;
LHSr(1:nflow,nflow+1:nshl*nflow) = Mr(1,2)*nr.*A1_r - Nr(1,2)*nr*K_r;
%% modification of LHS to account dirichlet BC
for i=1:nflow
    if(ibcg(1,i)== 1)
        LHSl(i,:) = 0.0;
    end
    if (ibcg(2,i)== 1)
        LHSr(i,:) = 0.0;
    end
end

end
