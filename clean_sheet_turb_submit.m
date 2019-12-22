% Heat equation model to understand Nusselt number in turbulent flow
% through sheet
% Can set up for heated walls or with internal dissipation profile

clear all
close all

Re_list=[10000];% 20000 30000 40000 50000 60000 70000 80000 90000 100000];
for Rei=1:length(Re_list)
    
    h=0.01;    % Subglacial gap half-width (m)
    
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
    y=0:dz:h;
    Re=Re_list(Rei); % Prescribe bulk Reynolds number
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
    
    % Dissipation profile (based on Abe and Antonia, missing interior portion)
    % depth-integrated E=<epsilon_bar>+<nu*(du/dz)^2>=(2.54.*log(yplus)+2.41).*u_tau^3;
    % eps here epsilon_bar, NOT depth-integrated!
    outer=find((h-abs(z))/h>0.2);
    eps(outer)=(2.45./((h-abs(z(outer)))./h)-1.7).*u_tau^3/h;
    inner=find((h-abs(z))/h<=0.2);
    eps(inner)=(2.54./((h-abs(z(inner)))./h)-2.6).*u_tau^3/h;
    % Inner portion of dissipation function (very near wall)
    hplus=u_tau*h/nu;
    [ii,jj]=find(abs(y/h-30/hplus)==(min(abs(y/h-30/hplus))));
    wall=find((h-abs(z))/h<30/hplus);
    eps(wall)=eps(jj);
    
    % % Get eddy diffusivity
    tau=tau_w.*z./h; % Laminar shear stress distribution - linear
    nu_T=-(tau(1:end-1)+tau(2:end))./2./(rhow.*diff(u')./dz)-nu; % Eddy viscosity (at half-nodes)
    nu_T(nu_T<0)=0;
    kappa_T=nu_T; % Eddy diffusivity from Reynolds' analogy (at half-nodes)
    
    dx=0.01;
    nx=100/dx;        % Number of x steps (iterations)
    
    alpha=(kappa+kappa_T').*dx/dz^2;
    
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
    
    if Q==0 % Heated walls
        % Type I at z=0
        a(1)=0;
        b(1)=1;
        c(1)=0;
        r(1)=1;
        
        
        % Type I at top
        a(J)=0;
        b(J)=1;
        c(J)=0;
        r(J)=1;
        
    elseif Q==1 % Dissipation case
        % Type I at z=0
        a(1)=0;
        b(1)=1;
        c(1)=0;
        r(1)=0;
        
        
        % Type I at top
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
            Q/rhow/cw*dx.*(eps(2:J-1)'.*rhow+mu.*((u(3:J)-u(1:J-2))./dz).^2);
        
        
        Tnext=thomas(a,b,c,r);
        T=Tnext';
        Tbar(it)=trapz(T)*dz/(2*h);
        
        % Save to plot for selected downstream distances
        if it*dx==1 | it*dx==5 | it*dx==10 | it*dx==20 | it*dx==50 | it*dx==100
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
    
    if Q==0 % For heated wall case (no dissipation)
        p=polyfit(x,log(T(1)-Tbar),1);
        Pe=trapz(u)*dz/(2*h)*4*h/(k/rhow/cw);
        Nu=-p(1)*Pe*(2*h)/2;
    elseif Q==1 % For internal dissipation case -- with constant T=0 walls
        Phi=(mu*trapz(((u(2:J)-u(1:J-1))./dz).^2)+trapz(eps(2:J-1))*rhow)*dz; % Integral of dissipation term - A/A profile with epsilon bar
        h_coeff=Phi/(2*(Tbar(end)-T(J)));
        Nu=h_coeff*4*h/k;
        E=(2.54*log(u_tau*h/nu)+2.41)*u_tau^3;  % This should hold for hplus>=300?  Compare E to Phi/rhow
        dissipation_integral=(trapz(eps(2:J-1))+trapz(((u(2:J)-u(1:J-1))./dz).^2)*nu)*dz/2; % <epsilon_bar> + <nu*(du/dy)^2> = u_tau^2*u_b (lhs integrates only over half-width, that's why /2)
        u_tau2u_b=u_tau^2*u_b;
    end
    
    % Comparison with Siegel and Sparrow for internal heat generation (fig. 2)
    SS=(T(1)-Tbar(end))./(Q*h^2/k); % Siegel and Sparrow, fig. 2
    
    % load AA_dissipationintegral_data
    % rhs=[rhs;u_tau2u_b];
    % lhs=[lhs;dissipation_integral];
    % Nu_list=[Nu_list;Nu];
    % save('AA_dissipationintegral_data','rhs','lhs','Nu_list')
    
    % load channel_heatedwalls_data
    % Nu_list=[Nu_list;Nu];
    % save('channel_heatedwalls_data','Nu_list')
    
end



