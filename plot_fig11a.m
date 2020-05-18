% Plot Figure 11
clear all;close all;

load Re1e5_jokhul_r2.mat

figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:nx/2),log(Tbar(1:nx/2)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all;load Re1e5_heatedwalls_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:nx/2),log(Tbar(1:nx/2)-T(end)),'--','linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all;load Re1e5_dissipation_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:nx/2),log(Tbar(1:nx/2)-T(end)),'--','linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

figure(2);legend('T_0=2 case','Heated wall case','Dissipation case','fontsize',14)
% title('Re=10^6','fontsize',16)