clear
clc
close all
%
global nflow
nflow = 3;
global Nnode
% geometry and mesh info
geo_input(1,1) = 0.0; %left end location
geo_input(1,2) = 0.4; %right end location
geo_input(1,3) = 0.2; %interface location
N_input(1,1) = 40; %total #of elements if no interface
N_input(1,2) = 20; %total #of elements in phase 0 if interface exit
N_input(1,3) = 20; %total #of elements in phase 1 if interface exit
ioption = 1; % DG formulation?1-ON, 2-OFF(single phase)
isolid = 0; % Solid formulation for single phase
            % 1-ON; 0-OFF
itau_f = 1; % stablization term for fluid
itau_s = -1; % stablization term for solid
%
NT = 100; %total time steps
CFL = 1;% CFL 
[Y_sol,x,mov] = FEM_1D(geo_input,N_input,ioption, NT, CFL, isolid,...
                       itau_f, itau_s);
% [Y_sol,x,mov] = FEM_1D_midpoint(geo_input,N_input,ioption, NT, CFL);
%
   if (ioption == 1)
     N0 = N_input(1,2);
     N1 = N_input(1,3);
     nplot0p = linspace(1,N0*nflow+1, N0+1);
     nplot0v = linspace(2,N0*nflow+2, N0+1);
     nplot0T = linspace(3,N0*nflow+3, N0+1);
%     
     nplot1p = linspace((N0+1)*nflow+1, (Nnode-1)*nflow+1,N1+1 );
     nplot1v = linspace((N0+1)*nflow+2, (Nnode-1)*nflow+2,N1+1 );
     nplot1T = linspace((N0+1)*nflow+3, (Nnode-1)*nflow+3,N1+1 );
%     
     figure
     plot(x(1:N0+1),Y_sol(nplot0v),'-+')
     hold on;
     pause(1)
     plot(x(N0+2:Nnode), Y_sol(nplot1v),'-*');
     legend ('air phase','air phase')
     pause(1)
     title ('Velocity')
%
     figure
     plot(x(1:N0+1),Y_sol(nplot0p),'-+')
     hold on;
     pause(1)
     plot(x(N0+2:Nnode), Y_sol(nplot1p),'-*');
     legend ('air phase','air phase')
     pause(1)
     title('Pressure')
%
     figure
     plot(x(1:N0+1),Y_sol(nplot0T),'-+')
     hold on;
     pause(1)
     plot(x(N0+2:Nnode), Y_sol(nplot1T),'-*');
     legend ('air phase','air phase')
     pause(1)
     title('Temperature')
   else
     N = N_input(1,1);
     nplotp = linspace(1,N*nflow+1, N+1);
     nplotv = linspace(2,N*nflow+2, N+1);
     nplotT = linspace(3,N*nflow+3, N+1);
     %
     figure
     plot(x,Y_sol(nplotp,1),'o-')
     pause(1)
     legend('air','Location','Northoutside')
     title('pressure')
     %
     figure
     plot(x,Y_sol(nplotv,1),'rx-')
     pause(1)
     legend('air','Location','Northoutside')
     title('velocity')
     %
     figure
     plot(x,Y_sol(nplotT,1),'kh-')
     pause(1)
     legend('air','Location','Northoutside')
     title('Temperature')
   end
   %%
%    figure
% %    set(gca,'Xlim',[0 0.4],'Ylim',[0 1.2]);
%    movie(mov);
   %%
%    v = VideoWriter ('1D_acoustic_DG_air_water_L_P_L_V_CFL0.5_stab.mp4','MPEG-4');
%    v.FrameRate = 10;
%    open(v);
%    writeVideo(v,mov(:));
%    close(v);
% %      
     