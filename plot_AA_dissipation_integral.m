% Plot Fig. 6
% Abe/Antonia
% <epsilon_bar> + <nu*(du/dy)^2> = u_tau^2*u_b

clear all;close all

% Channel flow with dissipation profile from Abe/Antonia
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000];

% With finer resolution dz=0.00001*h
load AA_dissipationintegral_data.mat
    
figure;plot(Re,lhs,Re,rhs,'o','linewidth',2)
xlabel('Re','fontsize',14)
set(gca,'fontsize',14)
legend('E=<{\Phi}>/2{\rho}','U_{\tau}^2U_b')