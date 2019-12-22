% Plot fully developed turbulent velocity profile (Figure 3)

Re_list=[10000 50000 100000];
for Rei=1:length(Re_list)

r0=0.01;    % Pipe radius (m)

% Typical 0 deg conditions used in ice sheet subglacial modeling
k=0.558; % Thermal conductivity of water (W/m/K)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
mu=1.787e-3; % Dynamic viscosity of water at 0 degC (Pa s)

rhow=1000; % Density of water (kg/m3)
kappa=k/(rhow*cw); % Thermal diffusivity
nu=mu/rhow; % Kinematic viscosity (m2/s)

Q=1; % Set to 0 to turn off internal heat generation, 1 for internal dissipation

dr=0.00001*r0;
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

figure(1);plot(u./u_tau,z./r0,'linewidth',2);hold on

end

set(gca,'fontsize',14)
xlabel('u/u_{\tau}','fontsize',14);ylabel('z/h or r/r_0','fontsize',14)
legend('Re=10^4','Re=5x10^4','Re=10^5')