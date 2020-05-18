% Plot Figure 11
clear all;close all;

load Re1e4_jokhul_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:15000),log(Tbar(1:15000)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all
load Re5e4_jokhul_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:15000),log(Tbar(1:15000)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all
load Re1e5_jokhul_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:15000),log(Tbar(1:15000)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all
load Re5e5_jokhul_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:15000),log(Tbar(1:15000)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

clear all
load Re1e6_jokhul_r2.mat
figure(2);hold on;
x=(0:(nx-1)).*dx;
plot(x(1:15000),log(Tbar(1:15000)-T(end)),'linewidth',2);
xlabel('x (m)','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)

figure(2);legend('Re=10^4','Re=5x10^4','Re=10^5','Re=5x10^5','Re=10^6','fontsize',14)
% title('Re=10^6','fontsize',16)