% Heat equation in radial coordinates (circular pipe) to understand Nusselt
% number in laminar flow
% Can set up for heated walls or with internal dissipation profile

clear all
close all

r0=0.1;    % Pipe radius (m)

% Typical 0 deg conditions used in ice sheet subglacial modeling
k=0.558; % Thermal conductivity of water (W/m/K)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
mu=1.787e-3; % Dynamic viscosity of water at 0 degC (Pa s)

rhow=1000; % Density of water (kg/m3)
kappa=k/(rhow*cw); % Thermal diffusivity
nu=mu/rhow; % Kinematic viscosity (m2/s)

Q=1; % Set to 0 to turn off internal heat generation, 1 for internal dissipation

dr=0.001*r0;    
J=ceil(r0/dr+1);
z=0:dr:r0;

% Get velocity profile
% For laminar flow
Pe=10;
u=3/2*Pe*k/rhow/cw/(2*r0).*(1-z'.^2./(r0)^2);
Re=u.*r0.*rhow./mu; % Local Reynolds' number
Re=mean(Re);

u=2.*1e-5.*(1-z'.^2./r0^2); % From textbook, test - circular laminar flow
ub=2/(r0^2)*trapz(u.*z')*dr; % Mean velocity

dx=0.001;
nx=3/dx;        % Number of x steps (iterations)

z_half=(z(1:end-1)+z(2:end))./2; % Half-node radial coordinates
alpha=(kappa).*dx./dr'.^2.*z_half';
% alpha=dx/dz^2/Pe*2/3.*ones(J,1);

% Define vectors (coefficients for interior nodes)
aconst=-alpha./2./z(2:end)';
bconst=(u(2:end-1)+alpha(2:end)./2./z(2:end-1)'+alpha(1:end-1)./2./z(2:end-1)');
cconst=-alpha./2./z(1:end-1)';

% Build a,b,c,r vectors
a=[0;aconst];
b=[0;bconst;0];
c=[cconst;0];
r=zeros(J,1);

% Boundary condition at pipe center (symmetry condition)
% Type II at z=0 (no-flux)
a(1)=0;
b(1)=-1;
c(1)=1;
r(1)=0;


if Q==0 % Type I at wall (heated)
    a(J)=0;
    b(J)=1;
    c(J)=0;
    r(J)=1;
end

if Q==1  % Type I at wall (T=0)
    a(J)=0;
    b(J)=1;
    c(J)=0;
    r(J)=0;
end

% % Type II at top (no-flux)
% a(J)=-1;
% b(J)=1;
% c(J)=0;
% r(J)=0;

% Initial condition u(x,0)=0
T=zeros(J,1);

T_CN=[];    % Pre-allocate matrix to store plot values for each t

% Loop through distance along x
for it=1:nx
    
    % Interior nodes
r(2:J-1)=alpha(1:end-1)./2./z(2:end-1)'.*T(1:J-2)+(u(2:J-1)-alpha(1:end-1)./2./z(2:end-1)'-alpha(2:end)./2./z(2:end-1)').*T(2:J-1)+alpha(2:end)./2./z(2:end-1)'.*T(3:J)+...
    Q/rhow/cw*dx.*(mu.*((u(3:J)-u(1:J-2))./dr).^2);


    Tnext=thomas(a,b,c,r);
    T=Tnext';
    Tbar(it)=2/(ub*r0^2)*trapz(u.*T.*z')*dr;
    
    % Save to plot for selected downstream distances
    if it*dx==0.01 | it*dx==0.1 | it*dx==0.5 | it*dx==1 | it*dx==2 | it*dx==3
    
        T_CN=[T_CN T];
    end
    
end


figure(1)
hold on
plot(T_CN,z./r0,'--','linewidth',1.2)
legend('x=0.01','x=0.1','x=0.2','x=0.5','x=1','x=2')
xlabel('T')
ylabel('r/r_{0}')

figure;
x=(0:(nx-1)).*dx;
plot(x,log(T(end)-Tbar));
xlabel('x');ylabel('log(T_w-T_{bar})')

% Nusselt number calculation
if Q==0 % For wall-heating case
    p=polyfit(x,log(T(end)-Tbar),1);
    Nu=-p(1)*rhow*cw*ub*r0^2/k;
elseif Q==1  % For internal dissipation case -- with constant T=0 walls
    Phi=(mu*2*pi*trapz(((u(2:J)-u(1:J-1))./dr).^2.*(z(2:J)+z(1:J-1))'/2))*dr;
    h_coeff=Phi/(2*pi*r0*(Tbar(end)-T(J)));
    Nu=h_coeff*2*r0/k;
end

Nu

% Comparison with Siegel and Sparrow for internal heat generation (fig. 2)
% -- make sure to set wall b.c. as no-flux for comparison
SS=(T(J)-Tbar(end))./(Q*r0^2/k); % Siegel and Sparrow, fig. 2



