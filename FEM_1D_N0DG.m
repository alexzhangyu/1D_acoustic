% % matlab code for 1-D acoustic problem
% % for research
function [Y_sol,x] = FEM_1D_N0DG( x0, x1,xif, N0,N1,NT,CFL )
% % x0: left end
% % x1: right end
% % xif: interface location
% % N0: number of elements in phase0
% % N1: number of elements in phase1
% % NT: total number of time steps
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % setting global parameters
global h0
global h1
%
h0 = ( xif - x0)/N0; % mesh size in phase0
h1 = ( x1 - xif)/N1; % mesh size in phase1
global dt
c_ss = 126;
c_air = 340;
dt = min(h0,h1)/max(c_air,c_air)*CFL;
% dt = 1.4e-5;
% total elements
global Nel
Nel = N0 + N1;
% total nodes
global Nnode
Nnode = N0 + N1 + 1; %accouting the duplicated node 
                     % at interface
% total number of equations
global nflow
nflow = 3;
% local degree of freedom
global nshl
nshl = 2;
% averaged main stream parameters
global p0
global T0
p0 = 1e5;
T0 = 288;
% time integation for flow
global am
global af
global gama
rinf = 0.7;
am = 0.5*(3.0 - rinf)/(1+rinf);
af = 1.0/(1+rinf);
gama = 0.5 + am -af;
% time integration for solid
global amb
global afb
global gamab
amb = 1.0;
afb = amb;
gamab = amb;
% % % % % material properties
% fluid %assuming air
global mu_l
global k_l
mu_l = 1.8e-5;
k_l = 0.026;
R = 8.314e3;
mv = 29; %molecular weight for air g/mol
rr = 1.4;
% solid
global bk_mod
global sh_mod
global E_mod
global k_s
bk_mod = 10.0e6;
sh_mod = 15.0e6;
E_mod = 2.0*sh_mod;
k_s = 31.1;
% common propterties: aa:related to temp
%                     bb: related to pressure
%                     rho_0 : intial density
%                     cv : specific heat per vol
% in the common properties arrays, 1st entry: fluid
%                                  2nd entry: solid
global rho_0
global aa
global bb
global cv_0
rho_0 = zeros(2,1);
aa = zeros(2,1);
bb = zeros(2,1);
cv_0 = zeros(2,1);
% 
rho_0(1,1) = 1.2111;
rho_0(2,1) = 1.8775e3;
% 
aa(1,1) = -1.0/T0;
aa(2,1) = -3.0*2.22e-5;
% 
bb(1,1) = 1.0/p0;
bb(2,1) = 1.0/bk_mod;
%
cv_0(1,1) = R/mv/(rr-1.0);
cv_0(2,1) = 1586;
% longitude wave speed for solid
 c_s = sqrt(E_mod/rho_0(2,1));
% %
% interface parameters
global S_k
global epln
global ita
global n0 
global n1
%% wave setup
lamda = 0.4;%wavelength
f = c_air/lamda;
global omg
omg = 2.0*pi*f;
%% some global arrays
% global coords
x = zeros(Nnode,1);
% filling global coords
for i = 1:N0+1
    x(i) = (i-1)*h0;
end
for i= N0+1:Nnode
    x(i) = xif + (i-(N0+1))*h1;
end
% global IEN array
IEN = zeros(Nel,3);
for i = 1:N0
    IEN(i,1) = i;
    IEN(i,2) = i+1;
    IEN(i,3) = 1; %solid phase
end
for i = (N0+1):Nel
    IEN(i,1)= i;
    IEN(i,2) = i+1;
    IEN(i,3) = 0; %fluid phase
end
% solution coeff vector Y, dot{Y} and
% Y_am, Y_af
Y =zeros(nflow*Nnode,1);
Y_dot = zeros(nflow*Nnode,1);
Y_am = zeros(nflow*Nnode,1);
Y_af = zeros(nflow*Nnode,1);
% solid strain vector
b11= zeros(Nel,1);
b11_dot =zeros(Nel,1);
% delta Y
sol = zeros(nflow*Nnode,1);
% GB
res = zeros(nflow*Nnode,1);
% global LHS
LHS = zeros(nflow*Nnode, nflow*Nnode);
%% apply the IC
% Y(:,1) = 1:nflow*Nnode; %%
Y(:,1) = 0.0;
Y_dot(:,1) = 0.0;
b11(:,1) = 0.0;
b11_dot(:,1) = 0.0;
Y_ini = Y;
Y_dot_ini =Y_dot;
%% time loop
for istp = 1:NT
%% enforing the BC
  Y = ensure_BC(Y,istp);
%flow predictor
  temp_Y = Y; %predicted y^{n+1}  
  temp_Y_dot = (temp_Y - Y)/(gama*dt)...
             +(1.0 - 1.0/gama)*Y_dot; % predicted y^dot^{n+1};
% flow corrector
   Y_af = Y + af*(temp_Y - Y);
   Y_am = Y_dot + am*(temp_Y_dot - Y_dot);
% enforcing the BC
   Y_af = ensure_BC(Y_af,istp);
% iteration loop
  for itr = 1:1
% loop over element   
   for iel = 1:Nel
%        Yl1 = zeros(nflow);
%        Yl2 = zeros(nflow);
% get the local Y from global Y,% could put in a loop and change the
% the data structure
       loc1s = (IEN(iel,1)-1)*nflow+1; % starting location for node1
       loc1e = IEN(iel,1)*nflow; % end location for node1
       loc2s = (IEN(iel,2)-1)*nflow+1; % starting location for node2
       loc2e = IEN(iel,2)*nflow; %end location for node2
% localized Y       
       Yaf_l1 = Y_af( loc1s : loc1e,1);
       Yaf_l2 = Y_af( loc2s : loc2e,1 );
       Yam_l1 = Y_am( loc1s : loc1e,1);
       Yam_l2 = Y_am( loc2s : loc2e,1);
        if(IEN(iel,3)==0)
            [resl,lhsl]=elm_mtx_liquid(Yaf_l1, Yaf_l2,...
                                       Yam_l1, Yam_l2,...
                                       aa(1,1),bb(1,1),...
                                       rho_0(1,1), cv_0(1,1),... 
                                       k_l, mu_l,...
                                       h1); %notice h1 here
        else
%              [resl,lhsl]=elm_mtx_solid(Yaf_l1, Yaf_l2,...
%                                        Yam_l1, Yam_l2,...
%                                        aa(2,1),bb(2,1),...
%                                        rho_0(2,1), cv_0(2,1),... 
%                                        k_s, sh_mod,...
%                                        h0,...
%                                        b11(iel,1),...
%                                        b11_dot(iel,1)); %notice h0 here
                [resl,lhsl]=elm_mtx_liquid(Yaf_l1, Yaf_l2,...
                                       Yam_l1, Yam_l2,...
                                       aa(1,1),bb(1,1),...
                                       rho_0(1,1), cv_0(1,1),... 
                                       k_l, mu_l,...
                                       h0); %notice h1 here

        end
 %% assembly to global array
 % assembly the residual % could improved by a loop
  res( loc1s : loc1e,1 ) = res( loc1s : loc1e,1)...
                           + resl(1:nflow,1);
  res( loc2s : loc2e,1 ) = res( loc2s : loc2e,1 )...
                           + resl(nflow+1:nshl*nflow,1);
 % assembly the global LHS % could be improved by a loop
   LHS(loc1s : loc1e, loc1s : loc1e) = LHS(loc1s : loc1e, loc1s : loc1e)...
                                     + lhsl( 1:nflow,1:nflow );
   LHS(loc1s : loc1e, loc2s : loc2e) = LHS(loc1s : loc1e, loc2s : loc2e)...
                                     + lhsl( 1:nflow,nflow+1:nshl*nflow );
   LHS(loc2s : loc2e, loc1s : loc1e) = LHS(loc2s : loc2e, loc1s : loc1e)...
                                     + lhsl( nflow+1:nshl*nflow,1:nflow );
   LHS(loc2s : loc2e, loc2s : loc2e) = LHS(loc2s : loc2e, loc2s : loc2e)...
                              + lhsl( nflow+1:nshl*nflow,nflow+1:nshl*nflow);
                                 
                                 
  
 
   end
%% applying the interface conditions
% some interface parameters
%
S_k = 1.0;
epln = 2.0e3;%need to be changed later
ita = zeros(3,3);
ita(2,2) = max(mu_l,E_mod/c_s);
ita(3,3) = max(k_l,k_s);
n0= 1.0; % normal
n1= -1.0; % normal
% get the local Y from global Y
locif0s = zeros(nshl); % size of nshl
locif0e = zeros(nshl);
locif1s = zeros(nshl);
locif1e = zeros(nshl);
%
for i = 1:nshl
    locif0s(i) = (IENIF0(1,i)-1)*nflow+1; % starting location for node0i
    locif0e(i) = IENIF0(1,i)*nflow; % end location for node0i
%    
    locif1s(i) = (IENIF1(1,i)-1)*nflow+1; % starting location for node1i
    locif1e(i) = IENIF1(1,i)*nflow; % end location for node1i
end
    Yaf_01 = Y_af( locif0s(1) : locif0e(1),1); % for 0 side
    Yaf_02 = Y_af( locif0s(2) : locif0e(2),1); % for 0 side
%    
    Yaf_11 = Y_af( locif1s(1) : locif1e(1),1); % for 1 side
    Yaf_12 = Y_af( locif1s(2) : locif1e(2),1); % for 1 side
%% resif0 = (nflow,1) !first nflow entries for 0 side or Gb_{2-}
%  resif1 = (nflow,1) ! next nflow entries for 1 side or Gb_{1+}
%  lhsif00 = (nflow, nflow)  !  d Gb_(2-)/ d Y_0
%  lhsif01 = (nflow, nflow)  !  d Gb_(2-)/ d Y_1
%  lhsif10 = (nflow, nflow)  !  d Gb_(1+)/ d Y_0
%  lhsif11 = (nflow, nflow)  !  d Gb_(1+)/ d Y_1
    [resif0, resif1, lhsif00, lhsif01, lhsif10, lhsif11] =...
                                 if_mtrx( Yaf_01, Yaf_02,...
                                          Yaf_11, Yaf_12,...
                                          b11(N0,1),...
                                          b11_dot(N0,1),...
                                          b11(N0+1,1),...
                                          b11_dot(N0+1,1) );
% assembly the interface conditons
   res(locif0s(2): locif0e(2)) = res(locif0s(2): locif0e(2))+...
                             resif0; % Gb_{2-}
   res(locif1s(1): locif1e(1)) = res(locif1s(1): locif1e(1))+...
                             resif1; % Gb_{1+}
% 
   for jclm = 1:nshl
     jtemp = nflow*(jclm-1)+1;  
     LHS(locif0s(2) : locif0e(2), locif0s(jclm) : locif0e(jclm)) =...
     LHS(locif0s(2) : locif0e(2), locif0s(jclm) : locif0e(jclm))...
                      + lhsif00( 1:nflow,jtemp:(jtemp+nflow-1) );
%                  
     LHS(locif0s(2) : locif0e(2), locif1s(jclm) : locif1e(jclm)) =...
     LHS(locif0s(2) : locif0e(2), locif1s(jclm) : locif1e(jclm))...
                      + lhsif01( 1:nflow,jtemp:(jtemp+nflow-1) );
%
     LHS(locif1s(1) : locif1e(1), locif1s(jclm) : locif1e(jclm)) =...
     LHS(locif1s(1) : locif1e(1), locif1s(jclm) : locif1e(jclm))...
                      + lhsif11( 1:nflow,jtemp:(jtemp+nflow-1) );
%
     LHS(locif1s(1) : locif1e(1), locif0s(jclm) : locif0e(jclm)) =...
     LHS(locif1s(1) : locif1e(1), locif0s(jclm) : locif0e(jclm))...
                      + lhsif10( 1:nflow,jtemp:(jtemp+nflow-1) );
     
   end
%% applying the boundary conditions (static)
   [LHS,res] = setBC(LHS,res);
%% linear solve by GMRES
%   sol = LHS\res;
    b = -res;
%     tol = 1e-9;
%     sol = gmres(LHS,b,100,tol);
%     b = -res(nflow+1:nflow*Nnode-2,1);
%     LHS_new = LHS(nflow+1:nflow*Nnode-2,nflow+1:nflow*Nnode-2);
% print out L2(b)
    norm (res,2)
    sol = LHS\b;
%     sol = LHS_new\b;
%% update Y_af
    Y_af = Y_af + sol;
    Y_am = (1.0 - am/gama)*Y_dot + am/dt/gama/af*(Y_af - Y);
%       Y_af(nflow+1:nflow*Nnode-2,1) = Y_af(nflow+1:nflow*Nnode-2,1) + sol;
%% enforcing the BC again 
    Y_af = ensure_BC(Y_af,istp);

  end %end of iteration  
    
%% update solid
    for iel = 1:N0 % hack, knowing where the solid is
       up1s = (IEN(iel,1)-1)*nflow+1; % starting location for node1
       up1e = IEN(iel,1)*nflow; % end location for node1
       up2s = (IEN(iel,2)-1)*nflow+1; % starting location for node2
       up2e = IEN(iel,2)*nflow; %end location for node2
% localized Y       
       Yaf_1s = Y_af( up1s : up1e,1);
       Yaf_2s = Y_af( up2s : up2e,1 );
% linearized approxi for grad Y
       gradv_af = (Yaf_2s(2,1)-Yaf_1s(2,1))/h0;
% update b11 and b11_dot
       b11(iel,1) = b11(iel,1) + dt*(amb-gamab)/amb*b11_dot(iel,1)...
                 + 2.0*gamab*dt/amb*gradv_af;
       b11_dot(iel,1) = (amb-1.0)/amb*b11_dot(iel,1)...
                    + 2.0/amb*gradv_af;
 
    end
%% update flow
     Y_old = Y;
     Y_dot_old = Y_dot;
     Y = Y_old + (Y_af - Y_old)/af;
%      Y_dot = (Y - Y_old)/(gama*dt)...
%              +(1.0 - 1.0/gama)*Y_dot_old;
       Y_dot = Y_dot_old + 1.0/am*(Y_am - Y_dot_old);

%% enforcing the BC again 
    Y = ensure_BC(Y,istp);
    Y_sol = Y;
%
     if (mod(istp,10) == 1)    
%% plot the solution
     ieq = 1;
     nplot0p = linspace(ieq,N0*nflow+ieq, N0+1);
     nplot0v = linspace(2,N0*nflow+2, N0+1);
     nplot0T = linspace(3,N0*nflow+3, N0+1);
%     
     nplot1p = linspace((N0+1)*nflow+ieq, (Nnode-1)*nflow+ieq,N1+1 );
     nplot1v = linspace((N0+1)*nflow+2, (Nnode-1)*nflow+2,N1+1 );
     nplot1T = linspace((N0+1)*nflow+3, (Nnode-1)*nflow+3,N1+1 );
%     
     figure
     plot(x(1:N0+1),Y(nplot0v),'-+')
     hold on;
     pause(1)
     plot(x(N0+2:Nnode), Y(nplot1v),'-*');
     legend ('air phase','air phase')
     pause(1)
     end  
  
 
end

%%


end

