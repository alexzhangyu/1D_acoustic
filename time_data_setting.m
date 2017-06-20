function time_data_setting(rinf)
% time integation for flow
global am
global af
global gama
am = 0.5*(3.0 - rinf)/(1+rinf);
af = 1.0/(1+rinf);
gama = 0.5 + am -af;
% time integration for solid
global amb
global afb
global gamab
% amb = 1.0;
% afb = amb;
% gamab = amb;
amb = am;
afb = af;
gamab = gama;
end