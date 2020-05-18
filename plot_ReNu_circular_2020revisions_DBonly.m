clear all;close all


% Circular pipe flow - wall heating - matches well with D-B correlation
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000];
Nu=[89.2587 165.3526 236.3955 304.3768 370.2091 434.4064 497.2896 559.0866 619.9480 680.0172 1.2511e+03 1.7895e+03 2.3086e+03 2.8093e+03 4.0257e+03 5.1997e+03]; % With finer mesh, from circular_heatedwalls_data
Pr=13.5146;
% figure(1);hold on;plot(Re,Nu,'-o','linewidth',2)
figure(1);semilogx(Re,Nu,'-o','linewidth',2)

% Circular pipe flow with internal dissipation with Abe/Antonia dissipation profile - epsilon*rhow in Phi
% calculation (compare Phi/rhow with E) - WITH NEAR-WALL EPSILON (CONSTANT from y/h<30/hplus)
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000];
Nu=[43.3028 85.510 126.5970 166.9143 206.6587 245.9381 284.8349 323.4069 361.6920 399.7326 763.6548 1.1172e+03 1.4616e+03 1.8287e+03 2.6788e+03 3.5119e+03]; % With finer mesh, from circular_dissipation_data
% figure(1);hold on;plot(Re,Nu,'-o','linewidth',2)
figure(1);hold on;semilogx(Re,Nu,'-o','linewidth',2)

Pr=13.5146;
Re2=10000:10000:1000000;
DB=0.024.*Re2.^0.8.*Pr^0.4;
Re2=10000:10000:1000000;
figure(1);hold on;semilogx(Re2,DB,'--','linewidth',2)

Gni=0.012.*(Re2.^0.87-280).*Pr^0.4;
figure(1);semilogx(Re2,Gni,'--','linewidth',2)

figure(1);set(gca,'FontSize',18)
axis([10000 1000000 0 6000])
xlabel('Re');ylabel('Nu');
legend('Heated Walls','Dissipation','Dittus-Boelter','Gnielinski')