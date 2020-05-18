% Heat equation in radial coordinates (circular pipe) to understand Nusselt number in turbulent flow
% Can set up for heated walls or with internal dissipation profile

clear all
close all

Re_list=[1000000];
Nu_list=zeros(size(Re_list));
for Rei=1:length(Re_list)
    
    r0=0.2;    % Pipe radius (m)
    
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
    
    ub=2/r0^2*trapz(u.*z')*dr; % Mean bulk velocity for circular pipe
    
    % Dissipation profile (based on Abe and Antonia, missing interior portion)
    % depth-integrated E=<epsilon_bar>+<nu*(du/dz)^2>=(2.54.*log(yplus)+2.41).*u_tau^3;
    % eps here epsilon_bar, NOT depth-integrated!
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
    
    dx=r0;
    nx=5000/dx;        % Number of x steps (iterations)
    
    z_half=(z(1:end-1)+z(2:end))./2; % Half-node radial coordinates
    alpha=(kappa+kappa_T').*dx./dr'.^2.*z_half';
    
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
    
    if Q==1 % Type I at wall (T=0)
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
    T=2*ones(J,1);
    
    T_CN=[];    % Pre-allocate matrix to store plot values for each t
    
    % Loop through distance along x
    for it=1:nx
        
        % Interior nodes
        r(2:J-1)=alpha(1:end-1)./2./z(2:end-1)'.*T(1:J-2)+(u(2:J-1)-alpha(1:end-1)./2./z(2:end-1)'-alpha(2:end)./2./z(2:end-1)').*T(2:J-1)+alpha(2:end)./2./z(2:end-1)'.*T(3:J)+...
            Q/rhow/cw*dx.*(eps(2:J-1)'.*rhow+mu.*((u(3:J)-u(1:J-2))./dr).^2);
        
        
        Tnext=thomas(a,b,c,r);
        T=Tnext';
        %     Tbar(it)=trapz(T)*dr/r0;
        Tbar(it)=2/(ub*r0^2)*trapz(u.*T.*z')*dr;
        
        % Save to plot for selected downstream distances
        if it*dx==1 | it*dx==5 | it*dx==10 | it*dx==20 | it*dx==50 | it*dx==100
            T_CN=[T_CN T];
        end
        
%         % Brinkman number
        Br_lam(it)=mu.*ub.^2./(k.*(T(end)-Tbar(it)));
        Phi=2*pi.*((mu*trapz(((u(2:J)-u(1:J-1))'./dr).^2.*(z(2:J)+z(1:J-1))/2))+trapz(eps.*z)*rhow)*dr;
        h_coeff=Phi/(2*pi*r0*(Tbar(end)-T(J)));
        Br(it)=Phi./(2*h_coeff*(2-T(J)));

        
    end
    
    % Plotting - turn off for timing in problem 3
    figure(1)
    hold on
    plot(T_CN,z./r0,'--','linewidth',1.2)
    % axis([0 1 0 1])
    legend('x=0.01','x=0.1','x=0.2','x=0.5','x=1','x=2','x=5','x=10','x=20','x=50','x=100','x=200','x=500','x=1000')
    xlabel('T')
    ylabel('r/r_{0}')
    
    figure(2);hold on;
    x=(0:(nx-1)).*dx;
    plot(x./r0,log(Tbar-T(end)),'linewidth',2);
    xlabel('x/r_0','fontsize',14);ylabel('ln(T_b-T_w)','fontsize',14)
    
    if Q==0 % For wall-heating case
        p=polyfit(x,log(T(end)-Tbar),1);
        Nu=-p(1)*rhow*cw*u_b*r0^2/k;
    elseif Q==1  % For internal dissipation case -- with constant T=0 walls
        Phi=2*pi.*((mu*trapz(((u(2:J)-u(1:J-1))'./dr).^2.*(z(2:J)+z(1:J-1))/2))+trapz(eps.*z)*rhow)*dr;
        h_coeff=Phi/(2*pi*r0*(Tbar(end)-T(J)));
        Nu=h_coeff*2*r0/k;
    end
    
    % Comparison with Siegel and Sparrow for internal heat generation (fig. 2)
    % -- make sure to set wall b.c. as no-flux for comparison
    SS=(T(J)-Tbar(end))./(Q*r0^2/k); % Siegel and Sparrow, fig. 2
    
    % Plot Brinkman number
    figure(3);hold on;plot(x./r0,Br)
    xlabel('x/r_0','fontsize',14);ylabel('Brinkman number')
    
    % Plot non-dimensional Brinkman
    figure(4);hold on;plot(x./r0,log((Tbar-T(J))./(2-T(J))),'linewidth',2)
    xlabel('x/r_0','fontsize',14);ylabel('ln((T_b-T_w)/(T_0-T_w))','fontsize',14)
    
end







