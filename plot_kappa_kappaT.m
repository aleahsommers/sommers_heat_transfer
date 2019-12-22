% Plot kappa/kappa_T (Figure 4)
clear all;close all

Re_list=[10000 50000 100000];
for Rei=1:length(Re_list)

r0=0.05;    % Pipe radius (m)

% Typical 0 deg conditions used in ice sheet subglacial modeling
k=0.558; % Thermal conductivity of water (W/m/K)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
mu=1.787e-3; % Dynamic viscosity of water at 0 degC (Pa s)

rhow=1000; % Density of water (kg/m3)
kappa=k/(rhow*cw); % Thermal diffusivity
nu=mu/rhow; % Kinematic viscosity (m2/s)

Q=1; % Set to 0 to turn off internal heat generation, 1 for internal dissipation

dr=0.01*r0;
J=ceil(r0/dr+1);
z=0:dr:r0;

% Get velocity profile
y=0:dr:r0;
Re=Re_list(Rei); % Prescribe bulk Reynolds number
vonK=0.4; % Von Karman constant
B=5; % Constant to use in log-law
Cf=fzero(@(Cf) -sqrt(2/Cf)+1/vonK*log(Re*sqrt(Cf)/2/sqrt(2))-1/vonK+B,0.009); % Friction factor (from log-law, A/A eqn 1.5)
u_b=Re*nu/(2*r0); % Bulk mean velocity -- A/A uses 2h, where h is half-width, but Re should depend on hydraulic radius?
tau_w=Cf*rhow*u_b^2/2; % Wall shear stress
u_tau=(tau_w/rhow)^0.5; % Friction velocity (inner velocity scale)
yplus=y.*u_tau./nu; % Non-dimensional y coordinate near wall

uplus(yplus<=20)=yplus(yplus<=20)-1.2533e-4*yplus(yplus<=20).^4+3.9196e-6.*yplus(yplus<=20).^5; % Velocity smooth profile (transition at yplus=20)
uplus(yplus>20)=1/vonK*log(yplus(yplus>20))+B; % Log-layer
u=flip(uplus)'.*(u_tau);

ub=2/r0^2*trapz(u.*z')*dr; % Mean bulk velocity for circular pipe

% Dissipation profile (based on Abe and Antonia, missing interior portion)
outer=find((r0-abs(z))/r0>0.2);
eps(outer)=(2.45./((r0-abs(z(outer)))./r0)-1.7).*u_tau^3/r0;
inner=find((r0-abs(z))/r0<=0.2);
eps(inner)=(2.54./((r0-abs(z(inner)))./r0)-2.6).*u_tau^3/r0;
% Inner portion of dissipation function (very near wall)
hplus=u_tau*r0/nu;
[ii,jj]=find(abs(y/r0-30/hplus)==(min(abs(y/r0-30/hplus))));
wall=find((r0-abs(z))/r0<30/hplus);
eps(wall)=eps(min(wall)-1);

% % Get eddy diffusivity
tau=tau_w.*z./r0; % Laminar shear stress distribution - linear
nu_T=-(tau(1:end-1)+tau(2:end))./2./(rhow.*diff(u')./dr)-nu; % Eddy viscosity (at half-nodes)
nu_T(nu_T<0)=0;
kappa_T=nu_T; % Eddy diffusivity from Reynolds' analogy (at half-nodes)

figure(1);plot([0 kappa_T]./kappa,z./r0,'linewidth',2);hold on
end

set(gca,'fontsize',14)
xlabel('{\kappa}_T/{\kappa}','fontsize',14);ylabel('z/h or r/r_0','fontsize',14)
legend('Re=10^4','Re=5x10^4','Re=10^5')