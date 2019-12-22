clear all;close all

load turb_circular_heatedwalls
Tnondim=(1-T_CN)/1;  % (Tw-Tb)/(Tw-T0)
figure;subplot(1,2,1)
plot(Tnondim(:,1:6),z./r0,'linewidth',2);hold on
xlabel('(T_w-T)/(T_w-T_0)','fontsize',14)
ylabel('r/r_0 or z/h','fontsize',14)
title('Heated Walls')

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

load turb_channel_heatedwalls
Thalf=T_CN((J+1)/2:end,:);
Tnondim=(1-Thalf)/1;  % (Tw-Tb)/(Tw-T0)
plot(Tnondim(:,1:6),z((J+1)/2:end)/h,'--','linewidth',2);

set(gca,'fontsize',14)

load turb_circular_dissipation
for i=1:6
    Tnondim(:,i)=T_CN(:,i)./T_CN(1,6);
end
subplot(1,2,2);hold on
plot(Tnondim(:,1:6),z./r0,'linewidth',2);
xlabel('T/T_{\infty}','fontsize',14)
ylabel('r/r_0 or z/h','fontsize',14)
title('Dissipation')

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

load turb_channel_dissipation
Thalf=T_CN((J+1)/2:end,:);
for i=1:6
    Tnondim(:,i)=Thalf(:,i)./Thalf(1,6);
end
plot(Tnondim(:,1:6),z((J+1)/2:end)/h,'--','linewidth',2);

set(gca,'fontsize',14)

legend('x/r_0=100','x/r_0=500','x/r_0=1{\times}10^3','x/r_0=2{\times}10^3','x/r_0=5{\times}10^3','x/r_0=10{\times}10^3',...
    'x/h=100','x/h=500','x/h=1{\times}10^3','x/h=2{\times}10^3','x/h=5{\times}10^3','x/h=10{\times}10^3',...
    'Location','south')
