% % matlab code for 1-D acoustic problem
% % for research
function [Y_sol,x,MOV] = FEM_1D(geo_input,N_input,ioption,NT,CFL, isolid,...
                                itau_f, itau_s)
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
   [ibcg,ibcb,bcg,bcb] = gen_BC_flag();
%    bcg(1,1) = 0.0;
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
global k_w
global mu_w
%
material_setting
%% setting dt
global dt
global c_air
global c_w
% c_air = 340;
c_air = 3.399978948629522e2;


if (ioption == 1)
%   dt = min(h(1,1),h(2,1))/max(c_s,c_air)*CFL;
%% NOTICE, ONLY FOR THE CASE WITH WATER + AIR + DG
  dt = min(h(1,1),h(2,1))/max(c_w,c_air)*CFL;
%% debugging
% c_air = 2.873492036849818e+02;
%   dt = h(1,1)/c_air*CFL;
%     dt = h(1,1)/1.0*CFL;  
elseif((ioption == 2)&&(isolid == 1))
  dt = h(1,1)/c_s*CFL;  
else
% testing
conv = 1.0;
% conv = 2.873492036849819e2;
%   dt = h(1,1)/c_air*CFL;
  dt = h(1,1)/conv*CFL;
end
% total number of equations
global nflow
% local degree of freedom
global nshl
nshl = 2;
% seting time integration
global am
global af
global gama
global amb
global afb
global gamab
time_data_setting(0.7);
global omg
  set_sinBC
% solution coeff vector Y, dot{Y} and
% Y_am, Y_af
Y =zeros(nflow*Nnode,1);
Y_dot = zeros(nflow*Nnode,1);
Y_am = zeros(nflow*Nnode,1);
Y_af = zeros(nflow*Nnode,1);
% solid strain vector
b11= zeros(Nel,1);
b11_dot =zeros(Nel,1);
%% apply the IC
Y(:,1) = 0.0;
Y_dot(:,1) = 0.0;
b11(:,1) = 0.0;
b11_dot(:,1) = 0.0;
% apply the gauss distribution IC
% Y = set_gauss_IC(N_input,ioption, x, Y);
Y_ini = Y; %storage
Y_dot_ini =Y_dot;
%% imposing the BC
  Y = impose_BC(Y,0);
%% time loop
for istp = 1:NT   
%  bcg(1,2)=1.0*sin(omg*istp*dt);
    
%flow predictor
  pred_Y = Y; %predicted y^{n+1}
  pred_Y_dot = (pred_Y - Y)/(gama*dt)...
             +(1.0 - 1.0/gama)*Y_dot; % predicted y^dot^{n+1};       
%imposing BC
  pred_Y = impose_BC(pred_Y,istp);
% flow corrector
   Y_af = Y + af*(pred_Y - Y);
   Y_am = Y_dot + am*(pred_Y_dot - Y_dot);
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
                                       h(1,1), itau_f); %notice h here
        else
%              [resl,lhsl]=elm_mtx_solid(Yaf_l1, Yaf_l2,...
%                                        Yam_l1, Yam_l2,...
%                                        aa(2,1),bb(2,1),...
%                                        rho_0(2,1), cv_0(2,1),... 
%                                        k_s, sh_mod,...
%                                        h(2,1),...
%                                        b11(iel,1),...
%                                        b11_dot(iel,1),...
%                                        itau_s); %notice h here
           [resl,lhsl]=elm_mtx_water(Yaf_l1, Yaf_l2,...
                                       Yam_l1, Yam_l2,...
                                       aa(3,1),bb(3,1),...
                                       rho_0(3,1), cv_0(3,1),... 
                                       k_w, mu_w,...
                                       h(2,1), itau_f); %notice h here 

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
%% modifing the LHS and/or RHS if no BC assigned at that point
   [resl,resr,LHSl,LHSr] = set_BCB(Y_af(1:3,1),Y_af(4:6,1),...
                                  Y_af(nflow*Nnode-5:nflow*Nnode-3,1),...
                                  Y_af(nflow*Nnode-2:nflow*Nnode,1),...
                                  IEN(1,3), IEN(Nel,3),...
                                  b11(1,1),  b11_dot(1,1),...
                                  b11(Nel,1),b11_dot(Nel,1),...
                                  h(1,1),    h(2,1));
    res(1:nflow,1) = res(1:nflow)+resl(1:nflow,1); % on left end
    res(nflow*Nnode-2:nflow*Nnode,1) = res(nflow*Nnode-2:nflow*Nnode,1)...
                                     + resr(1:nflow,1); % on right end
%
    LHS(1:nflow,1:nflow)= LHS(1:nflow,1:nflow) + LHSl(1:nflow, 1:nflow);
%
    LHS(1:nflow,nflow+1:nshl*nflow) = LHS(1:nflow,nflow+1:nshl*nflow)...
                                    + LHSl(1:nflow,nflow+1:nshl*nflow);
%
    LHS(nflow*Nnode-2:nflow*Nnode,nflow*Nnode-5:nflow*Nnode-3) =...
    LHS(nflow*Nnode-2:nflow*Nnode,nflow*Nnode-5:nflow*Nnode-3)...
    + LHSr(1:nflow, 1:nflow);
%
    LHS(nflow*Nnode-2:nflow*Nnode,nflow*Nnode-2:nflow*Nnode) =...
    LHS(nflow*Nnode-2:nflow*Nnode,nflow*Nnode-2:nflow*Nnode)...
    + LHSr(1:nflow, nflow+1:nshl*nflow);     
                                 
%% applying the boundary conditions
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
%     sol = linsolve(LHS,b);
      sol = LHS\b;
%     sol = LHS_new\b;
%% update Y_af
    Y_af = Y_af + sol;
    Y_am = Y_am + am/(dt*gama*af)*sol;
%     Y_am = (1.0 - am/gama)*Y_dot + am/(dt*gama*af)*(Y_af - Y);
%       Y_af(nflow+1:nflow*Nnode-2,1) = Y_af(nflow+1:nflow*Nnode-2,1) + sol;
%% enforcing the BC again 
%     Y_af = impose_BC(Y_af,istp);

  end %end of iteration  
    
%% update solid
  if (ioption == 1)
%     for iel = 1:N0 % hack, knowing where the solid is
%        up1s = (IEN(iel,1)-1)*nflow+1; % starting location for node1
%        up1e = IEN(iel,1)*nflow; % end location for node1
%        up2s = (IEN(iel,2)-1)*nflow+1; % starting location for node2
%        up2e = IEN(iel,2)*nflow; %end location for node2
% % localized Y       
%        Yaf_1s = Y_af( up1s : up1e,1);
%        Yaf_2s = Y_af( up2s : up2e,1 );
% % linearized approxi for grad Y
%        gradv_af = (Yaf_2s(2,1)-Yaf_1s(2,1))/h(2,1);
% % update b11 and b11_dot
%        b11(iel,1) = b11(iel,1) + dt*(amb-gamab)/amb*b11_dot(iel,1)...
%                  + 2.0*gamab*dt/amb*gradv_af;
%        b11_dot(iel,1) = (amb-1.0)/amb*b11_dot(iel,1)...
%                     + 2.0/amb*gradv_af;
%  
%     end
%%%%%
%       if (IEN(1,3) == 1)
%           [b11(1:N0,:),b11_dot(1:N0,:)] = update_solid( Y_af,N0,IEN(1:N0,:),...
%                                         b11(1:N0,:),  b11_dot(1:N0,:),...
%                                         h(1,1)); % default solid on the left
%       elseif (IEN(1,3) == 0)
%           [b11(N0+1:Nel,:),b11_dot(N0+1:Nel,:)] = update_solid( Y_af,N1,IEN(N0+1:Nel,:),...
%                                         b11(N0+1:Nel,:),  b11_dot(N0+1:Nel,:),...
%                                         h(2,1)); % notice h here
%       else
%           disp('error: wrong IEN for multiphase')
%       end
% % SPECIFIC FOR TEST PURE SOLID+DG      
%       if (IEN(N0+1,3) == 1)
%          [b11(N0+1:Nel,:),b11_dot(N0+1:Nel,:)] = update_solid( Y_af,N1,IEN(N0+1:Nel,:),...
%                                         b11(N0+1:Nel,:),  b11_dot(N0+1:Nel,:),...
%                                         h(2,1));
%       end
%      
  elseif((ioption == 2)&&(isolid == 1))
      [b11,b11_dot] = update_solid( Y_af,Nel,IEN,...
                                    b11,  b11_dot,...
                                    h(1,1)); % here h(1,1) == h(2,1)
  else
      disp('No solid phase')
    
  end
%% update flow
     Y_old = Y;
     Y_dot_old = Y_dot;
     Y = Y_old + (Y_af - Y_old)/af;
%      Y_dot = (Y - Y_old)/(gama*dt)...
%              +(1.0 - 1.0/gama)*Y_dot_old;
     Y_dot = Y_dot_old + 1.0/am*(Y_am - Y_dot_old);

%% enforcing the BC again 
%     Y = impose_BC(Y);
    Y_sol = Y;
%
    out_freq = 1;     
    if (mod(istp,out_freq) == 0 )
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
             x_solid_plot = linspace(1,Nel,Nel);
        % 
%              figure
             subplot(3,1,1)
             plot(x(1:N0+1),Y(nplot0p),'-+',x(N0+2:Nnode), Y(nplot1p),'-+')
%              axis([0 0.4 -5 50]);
             title ('Pressure')
             legend ('air','water')
             grid on;
             drawnow
%              hold on
%              hold on;
             pause(0.01)
             subplot(3,1,2)
             plot(x(1:N0+1),Y(nplot0v),'-+',x(N0+2:Nnode), Y(nplot1v),'-+')
%              axis([0 0.4 -0.02 0.12]);
             title('Velocity')
%              legend ('solid phase','air phase')
             grid on;
             drawnow
%
             subplot(3,1,3)
             plot(x(1:N0+1),Y(nplot0T),'-+',x(N0+2:Nnode), Y(nplot1T),'-+')
%              axis([0 0.4 -0.01 1]);
             title('Temperature')
             grid on;
%              legend ('solid phase','air phase')
%              axis([0 3.2 -100 2000]);
             drawnow
%              hold on
%              hold on;
             pause(0.01)
%              hold on
%              hold on;           
             
%                  subplot( 4,1, 4);
%                  plot(x_solid_plot,b11,'bo');
%                  axis([1 Nel -10e-2 10e-2]);
%                   title B11;
%                   grid on;
       else
%              if(isolid == 1)
%                  n_sub = 4;
%              else
%                  n_sub = 3;
%              end
             nplotp = linspace(1,N*nflow+1, N+1);
             nplotv = linspace(2,N*nflow+2, N+1);
             nplotT = linspace(3,N*nflow+3, N+1);
             x_solid = linspace(1,Nel,Nel);
             subplot(3,1,1);
             set(gcf, 'Position', [100, 100, 1600, 800])
             plot(x,Y(nplotp),'+-');
%              axis([0 3.2 -2000 2000]);
             title Pressure
             grid on;
             subplot(3,1,2);
             plot(x,Y(nplotv),'r+-');
%              axis([0 3.2 -5.0 5.0]);
%              legend('air','Location','Northoutside');
             title velocity;
             grid on;
%              set(gca,'Xlim',[0 0.4],'Ylim',[-0.2 1.2]);
             subplot(3,1,3);
             plot(x,Y(nplotT),'b*-');
%              axis([0 3.2 -8e-1 8e-1]);
%              legend('air');
             title Temperature;
             grid on;
%              if (n_sub > 3)
%              subplot( 4,1, 4);
%                  plot(x_solid,b11,'bo');
% %                  axis([1 Nel -10e-2 10e-2]);
%                   title B11;
%                   grid on;
%              end
             drawnow;
       end
             set(gcf, 'Position', [100, 100, 1600, 800]);
             MOV(floor(istp/out_freq))= getframe(gcf); 
    end
end
 
%%
end

