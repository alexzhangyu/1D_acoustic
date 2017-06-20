function set_IF_para
%% set interface conditon parameters
% interface parameters
global S_k
global epln
global ita
global n0 
global n1
% materal para
global mu_l
global mu_w
global E_mod
global sh_mod
global c_s
global k_l
global k_s
global k_w
% global parameters
global h
%
    S_k = 1.0;
    epln = 0.1;%need to be changed later
%     epln = 1.0;
    ita = zeros(3,3);
%     ita(2,2) = max(mu_l,sh_mod/c_s*h(2,1));
      ita(2,2) = max(mu_l,mu_w);
%      ita(2,2) = mu_l;
%      ita(2,2) = sh_mod/c_s*h(2,1);
%     ita(3,3) = max(k_l,k_s);
      ita(3,3) = max(k_l,k_w)*100;
%     ita(3,3) = k_l*100;
%     ita(3,3) = k_s;
    n0= 1.0; % normal
    n1= -1.0; % normal
end