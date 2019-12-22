% Plot turbulent and viscous dissipation profiles (Fig. 5)

clear all
close all

Relist=[10000 50000 100000];

for i=1:length(Relist)

h=0.01;    % Subglacial gap half-width (m)

% Typical 0 deg conditions used in ice sheet subglacial modeling
k=0.558; % Thermal conductivity of water (W/m/K)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
mu=1.787e-3; % Dynamic viscosity of water at 0 degC (Pa s)

rhow=1000; % Density of water (kg/m3)
kappa=k/(rhow*cw); % Thermal diffusivity
nu=mu/rhow; % Kinematic viscosity (m2/s)

Q=1; % Set to 0 to turn off internal heat generation

dz=0.001*h;    
J=ceil(2*h/dz+1);
z=-h:dz:h;

% Get velocity profile
y=0:dz:h;
Re=Relist(i); % Prescribe bulk Reynolds number
vonK=0.4; % Von Karman constant
B=5; % Constant to use in log-law
Cf=fzero(@(Cf) -sqrt(2/Cf)+1/vonK*log(Re*sqrt(Cf)/2/sqrt(2))-1/vonK+B,0.009); % Friction factor (from log-law, A/A eqn 1.5)
u_b=Re*nu/(2*h); % Bulk mean velocity -- A/A uses 2h, where h is half-width, but Re should depend on hydraulic radius?
tau_w=Cf*rhow*u_b^2/2; % Wall shear stress
u_tau=(tau_w/rhow)^0.5; % Friction velocity (inner velocity scale)
yplus=y.*u_tau./nu; % Non-dimensional y coordinate near wall
eta=y./h; % Non-dimensional y coordinate far from the wall
wake=2*pi/vonK*(3.*eta.^2-2.*eta.^3); % Wake function for defect layer
uplus(yplus<=20)=yplus(yplus<=20)-1.2533e-4*yplus(yplus<=20).^4+3.9196e-6.*yplus(yplus<=20).^5; % Velocity smooth profile (transition at yplus=20)
uplus(yplus>20)=1/vonK*log(yplus(yplus>20))+B; % Log-layer
u=[uplus'.*(u_tau);flip(uplus(1:end-1))'.*(u_tau)]; 

% Turbulent dissipation profile (based on Abe and Antonia, missing interior portion)
% depth-integrated E=<epsilon_bar>+<nu*(du/dz)^2>=(2.54.*log(yplus)+2.41).*u_tau^3;
% eps here epsilon_bar, NOT depth-integrated!
outer=find((h-abs(z))/h>0.2);
eps(outer)=(2.45./((h-abs(z(outer)))./h)-1.7).*u_tau^3/h;
inner=find((h-abs(z))/h<=0.2);
eps(inner)=(2.54./((h-abs(z(inner)))./h)-2.6).*u_tau^3/h;
% Inner portion of turbulent dissipation function (very near wall)
hplus=u_tau*h/nu;
[ii,jj]=find(abs(y/h-30/hplus)==(min(abs(y/h-30/hplus))));
wall=find((h-abs(z))/h<30/hplus);
eps(wall)=eps(jj);

vis=nu.*(diff(u)./dz).^2;
vis=[vis;vis(end)];

% % Get eddy diffusivity
tau=tau_w.*z./h; % Laminar shear stress distribution - linear
nu_T=-(tau(1:end-1)+tau(2:end))./2./(rhow.*diff(u')./dz)-nu; % Eddy viscosity (at half-nodes)
nu_T(nu_T<0)=0;
kappa_T=nu_T; % Eddy diffusivity from Reynolds' analogy (at half-nodes)


totaldiss=vis'+eps;

% Integrated total dissipation -- FOR CIRCULAR PIPE
r=0:dz:h;
epshalf=eps((J+1)/2:end); % Across radius r0 (center to wall)
vishalf=vis((J+1)/2:end)';
E=zeros(1,length(epshalf));epsint=zeros(1,length(epshalf));visint=zeros(1,length(epshalf));
E(1)=2*pi*r(1)*(epshalf(1)+vishalf(1));
epsint(1)=2*pi*r(1)*epshalf(1);
visint(1)=2*pi*r(1)*vishalf(1);
for i=2:length(epshalf)
    epsint(i)=trapz(epshalf(1:i))*2*pi*r(i)*dz;
    visint(i)=trapz(vishalf(1:i))*2*pi*r(i)*dz;
    E(i)=trapz(epshalf(1:i)+vishalf(1:i))*2*pi*r(i)*dz;
end
eps_norm=epsint./E;
vis_norm=visint./E;
E_norm=E./(E(end));


figure(1);subplot(1,2,1);plot(eps_norm,z((J+1)/2:end)/h,'linewidth',2);hold on;subplot(1,2,2);plot(vis_norm,z((J+1)/2:end)/h,'linewidth',2);hold on

end



figure(1);set(gca,'FontSize',18);
subplot(1,2,1);set(gca,'FontSize',18);xlabel('{\Phi}_T/{\Phi}');
ylabel('z/h or r/r_0');
 axis([0 1 0 1]);
subplot(1,2,2);xlabel('{\Phi}_{mean}/{\Phi}');axis([0 1 0 1]);
legend('Re=10^{4}','Re=5x10^{4}','Re=10^{5}')



