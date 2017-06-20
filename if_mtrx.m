function [resl0, resl1, lhsl00, lhsl01, lhsl10, lhsl11] = if_mtrx( yaf01, yaf02,...
                                                                   yaf11, yaf12,...
                                                                   b11_0,...
                                                                   b11_dot_0,...
                                                                   b11_1,...
                                                                   b11_dot_1,...
                                                                   IENIF0,...
                                                                   IENIF1,...
                                                                   h)
%
global nflow
global nshl
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
global cv_0
global aa
global bb
%
global S_k
global epln
global ita
global n0 
global n1
%
Nb02 = 1.0;
Nb11 = 1.0;
%
Nb02Na0 = [0.0 1.0];
Nb02Na1 = [1.0 0.0];
Nb02grad0 = [-0.5 0.5]*(2.0/h(1,1));
Nb02grad1 = [-0.5 0.5]*(2.0/h(2,1));
grad02Na0 = [0.0 0.5]*(2.0/h(1,1));
grad02Na1 = [0.5 0.0]*(2.0/h(1,1));
%
Nb11Na1 = [1.0 0.0];
Nb11Na0 = [0.0 1.0];
Nb11grad1 = [-0.5 0.5]*(2.0/h(2,1));
Nb11grad0 = [-0.5 0.5]*(2.0/h(1,1));
grad11Na1 = [-0.5 0.0]*(2.0/h(2,1));
grad11Na0 = [0.0 -0.5]*(2.0/h(2,1));
%
o02 = Nb02*n0;
o11 = Nb11*n1;
%setting all the sub matrix
c = [0.0 0.0 0.0;0.0 0.0 0.0;0.0 0.0 1.0];
%seting interface velocity v^{I}
%always using the phase 1 side velocity as v^{I}
% if(IENIF0(1,3) == 1)
%     vI = yaf02(2,1);
% elseif(IENIF0(1,3) == 0)
%     vI = yaf11(2,1);
% else
%     disp('Error: wrong material setting')
% end

% vI = 0.5*(yaf02(2,1) + yaf11(2,1));
%   vI = yaf02(2,1);
   vI = 0.0;
if(IENIF0(1,3) == 1)
  h0 = h(2,1);
% % A0
%   A0_0 = zeros(3,3);
%   A0_0(1,1) = rho_0(2,1)*bb(2,1);
%   A0_0(1,3) = rho_0(2,1)*aa(2,1);
%   A0_0(2,2) = rho_0(2,1);
%   A0_0(3,3) = rho_0(2,1)*cv_0(2,1);  
% %A1
%   A1_0 = zeros(3,3);
%   A1_0(1,2) = rho_0(2,1);
%   A1_0(2,1) = 1.0;
%   A1_0(3,2) = p0;
% %
%   K_0 = zeros(3,3);
%   K_0(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
%   K_0(3,3) = k_s;
% %
%   L_0 = zeros(3,1);
%   L_0(2,1) = (2.0/3.0)*sh_mod*( b11_0 + (dt*afb*(amb-gamab))/amb*b11_dot_0);
%% ONLY USED FOR THE TEST CASE WITH WATER+AIR+DG
  A1_0 = zeros(3,3);
  K_0 = zeros(3,3);
  L_0 = zeros(3,1);
  A0_0 = zeros(3,1);
%  
  [A0_0,A1_0,K_0] = set_elmmtx_water(A0_0,A1_0,K_0);
else
  h0 = h(1,1);
% A0
  A0_0 = zeros(3,3);
%   A0_0(1,1) = rho_0(1,1)*bb(1,1);
%   A0_0(1,3) = rho_0(1,1)*aa(1,1);
%   A0_0(2,2) = rho_0(1,1);
%   A0_0(3,3) = rho_0(1,1)*cv_0(1,1);
  A0_0 = setA0(A0_0, rho_0(1,1),aa(1,1),bb(1,1),cv_0(1,1) );
% A1
  A1_0 = zeros(3,3);
%   A1_0(1,2) = rho_0(1,1);
%   A1_0(2,1) = 1.0;
%   A1_0(3,2) = p0;
  A1_0 = setA1(A1_0, rho_0(1,1), p0);

%
  K_0 = zeros(3,3);
%   K_0(2,2) = (4.0/3.0)*mu_l;
%   K_0(3,3) = k_l;
  K_0 = setK_fluid(K_0, mu_l,k_l );
%
  L_0 = zeros(3,1);
end
%
if(IENIF1(1,3) == 1)
  h1 = h(2,1);
% % A0
%   A0_1 = zeros(3,3);
%   A0_1(1,1) = rho_0(2,1)*bb(2,1);
%   A0_1(1,3) = rho_0(2,1)*aa(2,1);
%   A0_1(2,2) = rho_0(2,1);
%   A0_1(3,3) = rho_0(2,1)*cv_0(2,1);    
% % A1
%   A1_1 = zeros(3,3);
%   A1_1(1,2) = rho_0(2,1);
%   A1_1(2,1) = 1.0;
%   A1_1(3,2) = p0;
% %
%   K_1 = zeros(3,3);
%   K_1(2,2) = (4.0/3.0)*sh_mod*(afb*gamab*dt)/amb;
%   K_1(3,3) = k_s;
% %
%   L_1 = zeros(3,1);
%   L_1(2,1) = (2.0/3.0)*sh_mod*( b11_1 + (dt*afb*(amb-gamab))/amb*b11_dot_1);
%% ONLY USED FOR THE TEST CASE WITH WATER+AIR+DG
  A1_1 = zeros(3,3);
  K_1 = zeros(3,3);
  L_1 = zeros(3,1);
  A0_1 = zeros(3,1);
%  
  [A0_1,A1_1,K_1] = set_elmmtx_water(A0_1,A1_1,K_1);

else
  h1 = h(1,1);
% A0
  A0_1 = zeros(3,3);
%   A0_1(1,1) = rho_0(1,1)*bb(1,1);
%   A0_1(1,3) = rho_0(1,1)*aa(1,1);
%   A0_1(2,2) = rho_0(1,1);
%   A0_1(3,3) = rho_0(1,1)*cv_0(1,1);
  A0_1 = setA0(A0_1, rho_0(1,1),aa(1,1),bb(1,1),cv_0(1,1) );  
% A1
  A1_1 = zeros(3,3);
%   A1_1(1,2) = rho_0(1,1);
%   A1_1(2,1) = 1.0;
%   A1_1(3,2) = p0;
  A1_1 = setA1(A1_1, rho_0(1,1), p0);
%
  K_1 = zeros(3,3);
%   K_1(2,2) = (4.0/3.0)*mu_l;
%   K_1(3,3) = k_l;
  K_1 = setK_fluid(K_1,mu_l,k_l);
%
  L_1 = zeros(3,1);
end
%
resl0 = zeros(nflow,1);
resl1 = zeros(nflow,1);
% local res
 resl0 = 0.5*Nb02Na0(1)*n0.*( A1_0*yaf01 )...
       + 0.5*Nb02Na0(2)*n0.*( A1_0*yaf02 );
 resl0 = resl0...
         -0.5*Nb02grad0(1)*n0.*(K_0*yaf01)...
         -0.5*Nb02grad0(2)*n0.*(K_0*yaf02);
 resl0 = resl0...
       + 0.5*S_k*grad02Na0(1)*n0.*(K_0*c*yaf01)...
       + 0.5*S_k*grad02Na0(2)*n0.*(K_0*c*yaf02);
 resl0 = resl0...
       + epln/h0*Nb02Na0(1)*n0*n0.*(ita*c*c*yaf01)...
       + epln/h0*Nb02Na0(2)*n0*n0.*(ita*c*c*yaf02);
% adding term 2 in the document   
 resl0 = resl0...
       + 0.5*Nb02Na0(1)*n0*vI.*(A0_0*yaf01)...
       + 0.5*Nb02Na0(2)*n0*vI.*(A0_0*yaf02);  
%   
  resl0 = resl0...
        + 0.5*Nb02Na1(1)*n0.*( A1_1*yaf11 )...
        + 0.5*Nb02Na1(2)*n0.*( A1_1*yaf12 );
  resl0 = resl0...
        - 0.5*Nb02grad1(1)*n0.*(K_1*yaf11)...
        - 0.5*Nb02grad1(2)*n0.*(K_1*yaf12);
  resl0 = resl0...
        + 0.5*S_k*grad02Na1(1)*n1.*(K_0*c*yaf11)...
        + 0.5*S_k*grad02Na1(2)*n1.*(K_0*c*yaf12);
  resl0 = resl0...
       + epln/h0*Nb02Na1(1)*n0*n1.*(ita*c*c*yaf11)...
       + epln/h0*Nb02Na1(2)*n0*n1.*(ita*c*c*yaf12);
% adding term 2 in the document
  resl0 = resl0...
       + 0.5*Nb02Na1(1)*n1*vI.*(A0_1*yaf11)...
       + 0.5*Nb02Na1(2)*n1*vI.*(A0_1*yaf12);
%   
  resl0 = resl0...
        - o02.* (L_0+L_1);
%%   
  resl1 = 0.5*Nb11Na1(1)*n1.*( A1_1*yaf11 )...
        + 0.5*Nb11Na1(2)*n1.*( A1_1*yaf12 );
  resl1 = resl1...
         -0.5*Nb11grad1(1)*n1.*(K_1*yaf11)...
         -0.5*Nb11grad1(2)*n1.*(K_1*yaf12);
  resl1 = resl1...
       + 0.5*S_k*grad11Na1(1)*n1.*(K_1*c*yaf11)...
       + 0.5*S_k*grad11Na1(2)*n1.*(K_1*c*yaf12);
  resl1 = resl1...
       + epln/h1*Nb11Na1(1)*n1*n1.*(ita*c*c*yaf11)...
       + epln/h1*Nb11Na1(2)*n1*n1.*(ita*c*c*yaf12);
% adding term 2 in the document
  resl1 = resl1...
       + 0.5*Nb11Na1(1)*n1*vI.*(A0_1*yaf11)...
       + 0.5*Nb11Na1(2)*n1*vI.*(A0_1*yaf12);
%
  resl1 = resl1...
        + 0.5*Nb11Na0(1)*n1.*( A1_0*yaf01 )...
        + 0.5*Nb11Na0(2)*n1.*( A1_0*yaf02 );
  resl1 = resl1...
        - 0.5*Nb11grad0(1)*n1.*(K_0*yaf01)...
        - 0.5*Nb11grad0(2)*n1.*(K_0*yaf02);
  resl1 = resl1...
       + 0.5*S_k*grad11Na0(1)*n0.*(K_1*c*yaf01)...
       + 0.5*S_k*grad11Na0(2)*n0.*(K_1*c*yaf02);
  resl1 = resl1...
       + epln/h1*Nb11Na0(1)*n1*n0.*(ita*c*c*yaf01)...
       + epln/h1*Nb11Na0(2)*n1*n0.*(ita*c*c*yaf02);
% adding term 2 in the document
  resl1 = resl1...
       + 0.5*Nb11Na0(1)*n0*vI.*(A0_0*yaf01)...
       + 0.5*Nb11Na0(2)*n0*vI.*(A0_0*yaf02);
%
  resl1 = resl1...
        - o11.* (L_0+L_1);
%% local lhs
   lhsl00 = zeros(nflow,nshl*nflow);
   lhsl01 = zeros(nflow,nshl*nflow);
   lhsl10 = zeros(nflow,nshl*nflow);
   lhsl11 = zeros(nflow,nshl*nflow);
%
   for jclm = 1: nshl
      jtemp = nflow*(jclm-1)+1;
      lhsl00(1:nflow,jtemp:(jtemp+nflow-1)) =...
          0.5*Nb02Na0(jclm)*n0.* A1_0...
         -0.5*Nb02grad0(jclm)*n0.* K_0...
         + 0.5*S_k*grad02Na0(jclm)*n0.*(K_0*c)...
         + epln/h0*Nb02Na0(jclm)*n0*n0.*(ita*c*c)...
         + 0.5*Nb02Na0(jclm)*n0*vI.*A0_0;
%
      lhsl01(1:nflow,jtemp:(jtemp+nflow-1)) =...
          0.5*Nb02Na1(jclm)*n0.*( A1_1)...
         -0.5*Nb02grad1(jclm)*n0.*(K_1)...
         + 0.5*S_k*grad02Na1(jclm)*n1.*(K_0*c)...
         + epln/h0*Nb02Na1(jclm)*n0*n1.*(ita*c*c)...
         + 0.5*Nb02Na1(jclm)*n1*vI.*A0_1;
%
      lhsl11(1:nflow,jtemp:(jtemp+nflow-1)) =...
          0.5*Nb11Na1(jclm)*n1.*( A1_1)...
         -0.5*Nb11grad1(jclm)*n1.*(K_1)...
         + 0.5*S_k*grad11Na1(jclm)*n1.*(K_1*c)...
         + epln/h1*Nb11Na1(jclm)*n1*n1.*(ita*c*c)...
         + 0.5*Nb11Na1(jclm)*n1*vI.* A0_1;
     
%
      lhsl10(1:nflow,jtemp:(jtemp+nflow-1)) =...
          0.5*Nb11Na0(jclm)*n1.*( A1_0 )...
         -0.5*Nb11grad0(jclm)*n1.*(K_0)...
         + 0.5*S_k*grad11Na0(jclm)*n0.*(K_1*c)...
         + epln/h1*Nb11Na0(jclm)*n1*n0.*(ita*c*c)...
         + 0.5*Nb11Na0(jclm)*n0*vI.* A0_0;
   end
   
        
         
end