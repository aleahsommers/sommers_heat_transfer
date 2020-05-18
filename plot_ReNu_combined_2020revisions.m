clear all;close all

% Circular pipe flow - wall heating - matches well with D-B correlation
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000];
Nu=[89.2587 165.3526 236.3955 304.3768 370.2091 434.4064 497.2896 559.0866 619.9480 680.0172 1.2511e+03 1.7895e+03 2.3086e+03 2.8093e+03 4.0257e+03 5.1997e+03]; % With finer mesh, from circular_heatedwalls_data
Pr=13.5146;
figure(1);semilogx(Re,Nu,'-o','linewidth',2)

% Circular pipe flow with internal dissipation with Abe/Antonia dissipation profile - epsilon*rhow in Phi
% calculation (compare Phi/rhow with E) - WITH NEAR-WALL EPSILON (CONSTANT from y/h<30/hplus)
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000];
Nu=[43.3028 85.510 126.5970 166.9143 206.6587 245.9381 284.8349 323.4069 361.6920 399.7326 763.6548 1.1172e+03 1.4616e+03 1.8287e+03 2.6788e+03 3.5119e+03]; % With finer mesh, from circular_dissipation_data
figure(1);hold on;semilogx(Re,Nu,'-o','linewidth',2)
    
% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

% Channel flow - wall heating
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000]; 
Nu=[162.4 304.0 436.5 563.6 686.8 807.1 925.0 1041 1155.3 1268.1 2.3430e+03 3.3585e+03 4.3391e+03 5.2859e+03 7.5885e+03 9.8136e+03]; % With finer mesh, from channel_heatedwalls_data

figure(1);hold on;semilogx(Re,Nu,'-.x','linewidth',2)
    
% Channel flow - with Abe/Antonia dissipation profile - epsilon*rhow in Phi
% calculation (compare Phi/rhow with E) - WITH NEAR-WALL EPSILON (CONSTANT from y/h<30/hplus)
Re=[10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 200000 300000 400000 500000 750000 1000000];
Nu=[91.1845 178.9022 264.1219 347.7000 430.0565 511.4583 592.0762 672.0347 751.4468 830.0000 1.5861e+03 2.3193e+03 3.0395e+03 3.8071e+03 5.5886e+03 7.3389e+03]; % With finer mesh, from AA_dissipationintegral_data
E=[0.0026 0.0173 0.0530 0.1176 0.2186 0.3632 0.5580 0.8099 1.1253 1.5105];
Phi=[4.98 33.29 101.75 224.15 415.74 684.34 1048 1512.7 2097.6 2794.3];
figure(1);semilogx(Re,Nu,'-.x','linewidth',2)
figure(1);set(gca,'FontSize',18)
axis([10000 1000000 0 10000])
xlabel('Re');ylabel('Nu');
legend('Circular conduit - heated walls','Circular conduit - dissipation','Sheet - heated walls','Sheet - dissipation')