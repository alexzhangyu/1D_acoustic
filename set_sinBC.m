function set_sinBC
global c_air
%% wave setup
lamda = 0.2;%wavelength
% f = c_air/lamda;
f = 0.5/lamda;
global omg
omg = 2.0*pi*f;

end
