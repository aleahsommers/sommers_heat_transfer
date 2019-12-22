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
    
% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

% Channel flow - wall heating
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000]; 
Nu=[162.4 304.0 436.5 563.6 686.8 807.1 925.0 1041 1155.3 1268.1]; % With finer mesh, from channel_heatedwalls_data

figure(1);hold on;plot(Re,Nu,'--o','linewidth',2)
    
% Channel flow - with Abe/Antonia dissipation profile - epsilon*rhow in Phi
% calculation (compare Phi/rhow with E) - WITH NEAR-WALL EPSILON (CONSTANT from y/h<30/hplus)
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000];
Nu=[91.1845 178.9022 264.1219 347.7000 430.0565 511.4583 592.0762 672.0347 751.4468 830.0000]; % With finer mesh, from AA_dissipationintegral_data
E=[0.0026 0.0173 0.0530 0.1176 0.2186 0.3632 0.5580 0.8099 1.1253 1.5105];
Phi=[4.98 33.29 101.75 224.15 415.74 684.34 1048 1512.7 2097.6 2794.3];
figure(1);plot(Re,Nu,'--o','linewidth',2)
figure(1);set(gca,'FontSize',18)
axis([10000 100000 0 1400])
xlabel('Re');ylabel('Nu');
legend('Circular conduit - heated walls','Circular conduit - dissipation','Sheet - heated walls','Sheet - dissipation')