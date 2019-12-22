% Heat equation laminar flow in sheet
% Aim to verify Dittus-Boelter or Gnielinski correlations for Nu
% Can set up for heated walls or with internal dissipation profile

clear all
close all

h=0.1;    % Subglacial gap half-width (m)

% Typical 0 deg conditions used in ice sheet subglacial modeling
k=0.558; % Thermal conductivity of water (W/m/K)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
mu=1.787e-3; % Dynamic viscosity of water at 0 degC (Pa s)

rhow=1000; % Density of water (kg/m3)
kappa=k/(rhow*cw); % Thermal diffusivity
nu=mu/rhow; % Kinematic viscosity (m2/s)

Q=1; % Set to 0 to turn off internal heat generation, 1 for internal dissipation

dz=0.001*h;
J=ceil(2*h/dz+1);
z=-h:dz:h;

% Get velocity profile
Pe=10;
u=3/2*Pe*k/rhow/cw/(2*h).*(1-z'.^2./(h)^2); % Channel laminar flow
Re=u.*h.*rhow./mu; % Local Reynolds' number
Re=mean(Re);

dx=0.001;
nx=3/dx;        % Number of x steps (iterations)

alpha=kappa.*dx/dz^2*ones(J-1,1);

% Define vectors (coefficients for interior nodes)
aconst=-alpha./2;
bconst=(u(2:end-1)+alpha(2:end)./2+alpha(1:end-1)./2);
cconst=-alpha./2;

% Build a,b,c,r vectors
a=[0;aconst];
b=[0;bconst;0];
c=[cconst;0];
r=zeros(J,1);

% Boundary conditions at walls
if Q==0 % Type I at walls (heated)
    a(1)=0;
    b(1)=1;
    c(1)=0;
    r(1)=1;
    
    a(J)=0;
    b(J)=1;
    c(J)=0;
    r(J)=1;
end

if Q==1 % Type I at walls (T=0)
    a(1)=0;
    b(1)=1;
    c(1)=0;
    r(1)=0;
    
    a(J)=0;
    b(J)=1;
    c(J)=0;
    r(J)=0;
end

% Initial condition u(x,0)=0
T=zeros(J,1);

T_CN=[];    % Pre-allocate matrix to store plot values for each t

% Loop through distance along x
for it=1:nx
    
    % Interior nodes
    r(2:J-1)=alpha(1:end-1)./2.*T(1:J-2)+(u(2:J-1)-alpha(1:end-1)./2-alpha(2:end)./2).*T(2:J-1)+alpha(2:end)./2.*T(3:J)+...
        Q/rhow/cw*dx.*(mu.*((u(3:J)-u(1:J-2))./dz).^2);
    
    
    Tnext=thomas(a,b,c,r);
    T=Tnext';
    Tbar(it)=trapz(T)*dz/(2*h);
    
    % Save to plot for selected downstream distances
    if it*dx==0.01 | it*dx==0.1 | it*dx==0.5 | it*dx==1 | it*dx==2 | it*dx==3
        
        T_CN=[T_CN T];
    end
    
end

figure(1)
hold on
plot(T_CN,z,'--','linewidth',1.2)
legend('x=0.01','x=0.1','x=0.2','x=0.5','x=1','x=2')
xlabel('T')
ylabel('z/b')

figure;
x=(0:(nx-1)).*dx;
plot(x,log(T(1)-Tbar));
xlabel('x');ylabel('log(T_w-T_{bar})')


if Q==1  % For internal dissipation case -- with constant T=0 walls
    Phi=(mu*trapz(((u(2:J)-u(1:J-1))./dz).^2))*dz;
    h_coeff=Phi/(2*(Tbar(end)-T(J)));
    Nu=h_coeff*4*h/k;
        
elseif Q==0  % For heated wall case (no dissipation)
    p=polyfit(x,log(T(1)-Tbar),1);
    Pe=trapz(u)*dz/(2*h)*4*h/(k/rhow/cw);
    Nu=-p(1)*Pe*(2*h)/2;
end

Nu


