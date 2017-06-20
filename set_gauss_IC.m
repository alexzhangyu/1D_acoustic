function [y] = set_gauss_IC( N_input,ioption, x, y )
% setting the IC of gaussian distribution

global nflow
global Nnode
% get the global dof which we would set the gauss IC on
if (ioption == 1)
    N0 = N_input(1,2);
    N1 = N_input(1,3); 
    ndof0p = linspace(1,N0*nflow+1, N0+1);
    ndof0v = linspace(2,N0*nflow+2, N0+1);
    ndof0T = linspace(3,N0*nflow+3, N0+1);
    %
    ndof1p = linspace((N0+1)*nflow+1, (Nnode-1)*nflow+1,N1+1 );
    ndof1v = linspace((N0+1)*nflow+2, (Nnode-1)*nflow+2,N1+1 );
    ndof1T = linspace((N0+1)*nflow+3, (Nnode-1)*nflow+3,N1+1 );
%
    dofloc0 = ndof0v;
    dofloc1 = ndof1v;
    x_endloc0 = N0+1;
    x_startloc1 = N0+2;
    mean = 1.0;
    dev = 0.03;
    y(dofloc0,1) = normpdf(x(1:x_endloc0,1),mean,dev);
%       y(dofloc0,1) = 10*sin(2*pi/1.6*x(1:x_endloc0,1));
%     y(dofloc1,1) = normpdf(x(x_startloc1:Nnode,1),mean,dev);
%     y(dofloc0,1) = 12.0;
else
    N = N_input(1,1);
    ndofp = linspace(1,N*nflow+1, N+1);
    ndofv = linspace(2,N*nflow+2, N+1);
    ndofT = linspace(3,N*nflow+3, N+1);
%
    dofloc = ndofv;
    x_endloc = Nnode;
    mean = 1.0;
    dev = 0.04;
    y(dofloc,1) = normpdf(x(1:x_endloc,1),mean,dev);
%     N0 = N_input(1,2);
%     N1 = N_input(1,3); 
%     ndof0p = linspace(1,N0*nflow+1, N0+1);
%     ndof0v = linspace(2,N0*nflow+2, N0+1);
%     ndof0T = linspace(3,N0*nflow+3, N0+1);
%     %
%     ndof1p = linspace((N0+1)*nflow+1, (Nnode-1)*nflow+1,N1+1 );
%     ndof1v = linspace((N0+1)*nflow+2, (Nnode-1)*nflow+2,N1+1 );
%     ndof1T = linspace((N0+1)*nflow+3, (Nnode-1)*nflow+3,N1+1 );
% %
%     dofloc0 = ndof0p;
%     dofloc1 = ndof1v;
%     x_endloc0 = N0+1;
%     x_startloc1 = N0+2;
%     mean = 1.2;
%     dev = 0.04;
%     y(dofloc0,1) = normpdf(x(1:x_endloc0,1),mean,dev);
end
   
end

