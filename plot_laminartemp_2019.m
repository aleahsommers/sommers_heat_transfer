clear all;close all

load laminar_circular_heatedwalls
Tnondim=(1-T_CN)/1;  % (Tw-Tb)/(Tw-T0)
figure;
subplot(1,2,1);
plot(Tnondim(:,2:5),z./r0,'linewidth',2);hold on
xlabel('(T_w-T)/(T_w-T_0)')
ylabel('r/r_0 or z/h')
title('Heated Walls')

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

load laminar_channel_heatedwalls
Thalf=T_CN((J+1)/2:end,:);
Tnondim=(1-Thalf)/1;  % (Tw-Tb)/(Tw-T0)
plot(Tnondim(:,2:5),z((J+1)/2:end)/h,'--','linewidth',2);
set(gca,'fontsize',14)

load laminar_circular_dissipation
for i=1:6
    Tnondim(:,i)=T_CN(:,i)./T_CN(1,6);
end
subplot(1,2,2);
plot(Tnondim(:,2:5),z./r0,'linewidth',2);hold on
axis([-0.5 1 0 1])
xlabel('T/T_{\infty}')
ylabel('r/r_0 or z/h')
title('Dissipation')

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

load laminar_channel_dissipation
Thalf=T_CN((J+1)/2:end,:);
for i=1:6
    Tnondim(:,i)=Thalf(:,i)./Thalf(1,6);
end
plot(Tnondim(:,2:5),z((J+1)/2:end)/h,'--','linewidth',2);
axis([0 1 0 1])
set(gca,'fontsize',14)
legend('x/r_0=1','x/r_0=5','x/r_0=10','x/r_0=20',...
    'x/h=1','x/h=5','x/h=10','x/h=20',...
    'Location','southeast')


