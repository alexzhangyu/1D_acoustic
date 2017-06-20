% % matlab code for 1-D acoustic problem
% % for research
function [Y_sol,x,MOV] = FEM_1D_midpoint(geo_input,N_input,ioption,NT,CFL )
% % ioption : 1-DG, 2-CG
% % NT: total number of time steps
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % setting geometry info
global h
global Nel
global Nnode

if (ioption == 1)
    x0 = geo_input(1,1);
    x1 = geo_input(1,2);
    xif = geo_input(1,3);
    N0 = N_input(1,2);
    N1 = N_input(1,3);
%       x0 = 0;
%     x1 = 0.4;
%     xif = 0.2;
%     N0 = 20;
%     N1 = 20;
    [IEN,IENIF0,IENIF1,x,h,Nel,Nnode] = geo_IF(x0, x1,xif, N0,N1);
else
    x0 = geo_input(1,1);
    x1 = geo_input(1,2);
    N = N_input(1,1);
    [IEN,x,h,Nel,Nnode] = geo_CG(x0, x1, N);
end
%% imposing BC flag and BC conditions
global ibcg
global bcg
global ibcb
global bcb
   [ibcg,ibcb,bcg,bcb] = gen_BC_flag;
%% setting mainstream
   set_p0T0();
%% material properties
%air
global mu_l
global k_l
% solid
global bk_mod
global sh_mod
global k_s
%shared
global rho_0
global aa
global bb
global cv_0
global c_s
%
material_setting
%% setting dt
global dt
global c_air
global conv
c_air = 340;

if (ioption == 1)
  dt = min(h(1,1),h(2,1))/max(c_s,c_air)*CFL;
else
% testing
conv = 1;
%   dt = h(1,1)/c_air*CFL;
  dt = h(1,1)/conv*CFL;
end
% total number of equations
global nflow
% local degree of freedom
global nshl
nshl = 2;
% % seting time integration
% global am
% global af
% global gama
% global amb
% global afb
% global gamab
% time_data_setting(0.6);
% set sin BC
  set_sinBC
% solution coeff vector Y, dot{Y} and
% Y_am, Y_af
Y =zeros(nflow*Nnode,1);
% Y_dot = zeros(nflow*Nnode,1);
% Y_am = zeros(nflow*Nnode,1);
% Y_af = zeros(nflow*Nnode,1);
% solid strain vector
b11= zeros(Nel,1);
b11_dot =zeros(Nel,1);
%% apply the IC
% Y(:,1) = 1:nflow*Nnode; %%
Y(:,1) = 0.0;
Y_dot(:,1) = 0.0;
b11(:,1) = 0.0;
b11_dot(:,1) = 0.0;
Y_ini = Y; %storage
Y_dot_ini =Y_dot;
%% imposing the BC
%   Y = impose_BC_midpoint(Y,0);
%% time loop
for istp = 1:NT     
% % enforcing the BC
%    Y = impose_BC_midpoint(Y,istp);
% storing Y_old
   Y_old = Y;
% iteration loop
  for itr = 1:1
% zero delta Y
sol = zeros(nflow*Nnode,1);
% zero GB
res = zeros(nflow*Nnode,1);
% zero global LHS
LHS = zeros(nflow*Nnode, nflow*Nnode);      
% loop over element   
   for iel = 1:Nel
% get the local Y from global Y,% could put in a loop and change the
% the data structure
       loc1s = (IEN(iel,1)-1)*nflow+1; % starting location for node1
       loc1e = IEN(iel,1)*nflow; % end location for node1
       loc2s = (IEN(iel,2)-1)*nflow+1; % starting location for node2
       loc2e = IEN(iel,2)*nflow; %end location for node2
% localized Y       
       Y_l1 = Y_old( loc1s : loc1e,1);
       Y_l2 = Y_old( loc2s : loc2e,1 );
%        Yam_l1 = Y_am( loc1s : loc1e,1);
%        Yam_l2 = Y_am( loc2s : loc2e,1);
        if(IEN(iel,3)==0)
            [resl,lhsl]=elm_mtx_midpoint(Y_l1, Y_l2,...
                                       aa(1,1),bb(1,1),...
                                       rho_0(1,1), cv_0(1,1),... 
                                       k_l, mu_l,...
                                       h(1,1)); %notice h here
        else
            [resl,lhsl]=elm_mtx_midpoint(Y_l1, Y_l2,...
                                       aa(1,1),bb(1,1),...
                                       rho_0(1,1), cv_0(1,1),... 
                                       k_l, mu_l,...
                                       h(1,1)); %notice h here


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
   if (ioption ==1)
%% applying the interface conditions
% setting interface parameters
    set_IF_para();
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
        b11_dot(N0+1,1),....
        IENIF0,...
        IENIF1,...
        h);
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
   end
%% applying the boundary conditions
   res = impose_BC_midpoint(res,istp);
   LHS = setLHS_midpoint(LHS);
%% linear solve by GMRES
%   sol = LHS\res;
    b = res;
%     tol = 1e-9;
%     sol = gmres(LHS,b,100,tol);
%     b = -res(nflow+1:nflow*Nnode-2,1);
%     LHS_new = LHS(nflow+1:nflow*Nnode-2,nflow+1:nflow*Nnode-2);
% print out L2(b)
    norm (res,2)
%     sol = linsolve(LHS,b);
      sol = LHS\b;
      Y = sol;
% debug
      Y = impose_BC_midpoint(Y,istp);
%     sol = LHS_new\b;

  end %end of iteration  
    
%% update solid
  if (ioption == 1)
    for iel = 1:N0 % hack, knowing where the solid is
       up1s = (IEN(iel,1)-1)*nflow+1; % starting location for node1
       up1e = IEN(iel,1)*nflow; % end location for node1
       up2s = (IEN(iel,2)-1)*nflow+1; % starting location for node2
       up2e = IEN(iel,2)*nflow; %end location for node2
% localized Y       
       Yaf_1s = Y_af( up1s : up1e,1);
       Yaf_2s = Y_af( up2s : up2e,1 );
% linearized approxi for grad Y
       gradv_af = (Yaf_2s(2,1)-Yaf_1s(2,1))/h(2,1);
% update b11 and b11_dot
       b11(iel,1) = b11(iel,1) + dt*(amb-gamab)/amb*b11_dot(iel,1)...
                 + 2.0*gamab*dt/amb*gradv_af;
       b11_dot(iel,1) = (amb-1.0)/amb*b11_dot(iel,1)...
                    + 2.0/amb*gradv_af;
 
    end
  end
%  
    Y_sol = Y;
%
      if (mod(istp,10) == 0)    
%% plot the solution
       if (ioption == 1)
             nplot0p = linspace(1,N0*nflow+1, N0+1);
             nplot0v = linspace(2,N0*nflow+2, N0+1);
             nplot0T = linspace(3,N0*nflow+3, N0+1);
        %     
             nplot1p = linspace((N0+1)*nflow+1, (Nnode-1)*nflow+1,N1+1 );
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
       else
             nplotp = linspace(1,N*nflow+1, N+1);
             nplotv = linspace(2,N*nflow+2, N+1);
             nplotT = linspace(3,N*nflow+3, N+1);
             plot(x,Y(nplotp),'o-');
%              axis([0 0.4 0 1.2]);
             legend('air');
             title pressure;
             grid on;
             set(gca,'Xlim',[0 0.4],'Ylim',[-0.2 1.2]);
             drawnow
       end
             MOV(istp/10)= getframe(gcf);
      end  
  
 
end
 
%%
end

