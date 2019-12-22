clear all;close all


% Circular pipe flow - wall heating - matches well with D-B correlation
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000];
Nu=[89.2587 165.3526 236.3955 304.3768 370.2091 434.4064 497.2896 559.0866 619.9480 680.0172]; % With finer mesh, from circular_heatedwalls_data
Pr=13.5146;
figure(1);hold on;plot(Re,Nu,'-o','linewidth',2)

% Circular pipe flow with internal dissipation with Abe/Antonia dissipation profile - epsilon*rhow in Phi
% calculation (compare Phi/rhow with E) - WITH NEAR-WALL EPSILON (CONSTANT from y/h<30/hplus)
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000];
Nu=[43.3028 85.510 126.5970 166.9143 206.6587 245.9381 284.8349 323.4069 361.6920 399.7326]; % With finer mesh, from circular_dissipation_data
figure(1);hold on;plot(Re,Nu,'-o','linewidth',2)

Pr=13.5146;
DB=0.024.*Re.^0.8.*Pr^0.4;
figure(1);hold on;plot(Re,DB,'linewidth',2)

figure(1);set(gca,'FontSize',18)
axis([10000 100000 0 800])
xlabel('Re');ylabel('Nu');
legend('Heated Walls','Dissipation','Dittus-Boelter')