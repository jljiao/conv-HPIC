%% Hybrid Particle-in-Cell (HPIC) Code （1D）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The HPIC1D is an one-dimentional hybrid PIC code, in which the ions are 
% treated kinetically and the electrons are assumed to be a massless fluid.
% This code is used to test a convolution HPIC method for non-quasineutral
% collisionless plasma electrostatic field.
% This code is written by J. L. Jiao.
clc;clear;
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In c.g.s Unit
e = 4.8e-10; % Elementary Charge.
mp = 1.67e-24; % Proton Mass [g].
c = 3e10; % Speed of Light [cm/s].
kB  = 1.38e-16; % Boltzmann Constant.
%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
t = 5e-12; % Total Simulation Time [s]
dt_factor = 0.1; % Ratio of dt and 1/omega_pi.
Lx = 1e-2; % The Length of the Simulation Box [cm].
NxDebye = 4; % The Cell Numbers in Debye Length.
Np = 1e6; % The Total Number of the Ion Pesudo-Particles. 
ND = 4; % ND*lambda_De is the Kernel Size

% Plasma Parameters
Z = 1; % Atom Number.
A = 1; % Mass Number.
ne0 = 1.0e20; % Electron Density [/cm^3].
Te0 = 40*1e3; % Electron Temperature [eV].
Ti0 = 1e3; % Ion Temperature [eV].
ux0 = 0.01*c; % Fluid Velocity in X Direction [cm/s].
gamma = 3; % gamma=(f+2)/f, f is free degree. For 1D case, f=1.

% Physics Option
phys = 2; % 1: HPIC.
          % 2: Convolution HPIC.
%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qi = Z*e; % Charge of Ion.
mi = mp*A; % Mass of Ion.
ni0 = ne0/Z; % Density of Ion.
lambda_De = 743*(Te0/ne0)^0.5; % The Debye Length [cm].
omega_pi = 1.32e3*Z*(ni0/A)^0.5; % Ion Plasma Frequency [/s].

dx = lambda_De/NxDebye; % Cell Size.
Nx = floor(Lx/dx); % Cell Numbers.
if(mod(Nx,2)==1)
    Nx = Nx+1;
end
Lx = Nx*dx;
ppc = floor(Np/Nx); % The Number of the Ion Pesudo-Particles per Cell. 
dt = 1/omega_pi*dt_factor; % Time Step Size.
Nt = floor(t/dt); % Numbers of Time Step.
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids Initialization
x = linspace(-1.5*dx,Lx+1.5*dx,Nx+4); % Coordinates of the Box.
ni = zeros(1,Nx+3); % Ion Density.
ne_prev = ones(1,Nx+3)*ne0; % Electron Density.
E = zeros(1,Nx+4); % Electrostatic Field.

% Particles Initialization
np = (Nx/2+1)*ppc;
xp = [rand(np,1)*(Lx/2+dx)-dx;...
    rand(np,1)*(Lx/2+dx)+Lx/2]; % X Position of the Ion Pesudo-Particles.
vxp = [randn(np,1)*(2*Ti0/A/9.38e8)^0.5*c+ux0;...
    randn(np,1)*(2*Ti0/A/9.38e8)^0.5*c-ux0]; % vx of the Ion Pesudo-Particles
wp = ni0/ppc*ones(2*np,1); % Weight of the Ion Pesudo-Particles
%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nt
    % Calculate Density
    ind0 = floor(xp/dx);
    ind1 = ind0+1;
    w0 = 1-(xp/dx-ind0);
    w1 = 1-w0;
    n0 = w0.*wp;
    n1 = w1.*wp;
    for j=1:length(xp)
        ni(ind0(j)+2) = ni(ind0(j)+2)+n0(j);
        ni(ind1(j)+2) = ni(ind1(j)+2)+n1(j);
    end
    ni(1) = ni0;
    ni(end) = ni0;
    ni = movmean(ni,3); % Smooth Density
    ne = Z*ni;
    
    % Calculate Electrostatic Field
    switch  phys
        case 1
            Te = Te0*(ne/ne0).^(gamma-1);
            lDe = 743*(Te./ne).^0.5;
            E(2:end-1) = -4*pi*e*gamma*((lDe(2:end)+lDe(1:end-1))/2).^2.*...
                (ne(2:end)-ne(1:end-1))/dx;
            E(1) = 0.0;
            E(end) = 0.0;
            E = movmean(E,3); % Smooth Electrostatic Field
        case 2
            Te = Te0*(ne_prev/ne0).^(gamma-1);
            lDe = 743*(Te./ne_prev).^0.5;
            ne_temp = [ne0*ones(1,ND*NxDebye),ne,ne0*ones(1,ND*NxDebye)];
            for j=1:Nx+3
                kernel = exp(-abs((-ND*NxDebye:1:ND*NxDebye)*dx/lDe(j)/gamma^0.5));
                kernel = kernel/sum(kernel);
                ne(j) = sum(ne_temp(j:j+2*ND*NxDebye).*kernel);
            end

            E(2:end-1) = -4*pi*e*gamma*((lDe(2:end)+lDe(1:end-1))/2).^2.*...
                (ne(2:end)-ne(1:end-1))/dx;
            E(1) = 0.0;
            E(end) = 0.0;
    end

    % Push Ions
    vxp = vxp+qi*dt/mi*interp1(x,E,xp);
    xp = xp+vxp*dt;
    
    % Injector Boundary for Right Side
    vxp(xp>=Lx) = [];
    wp(xp>=Lx) = [];
    xp(xp>=Lx) = [];
    
    xpbr = Lx+rand(ppc,1)*dx;
    vxpbr = randn(ppc,1)*(2*Ti0/A/9.38e8)^0.5*c-ux0;
    wpbr = ni0/ppc*ones(ppc,1);
    
    xp = [xp;xpbr];
    vxp = [vxp;vxpbr];
    wp = [wp;wpbr];
    
    % Injector Boundary for Left Side
    vxp(xp<=0) = [];
    wp(xp<=0) = [];
    xp(xp<=0) = [];
    
    xpbl = -rand(ppc,1)*dx;
    vxpbl = randn(ppc,1)*(2*Ti0/A/9.38e8)^0.5*c+ux0;
    wpbl = ni0/ppc*ones(ppc,1);
    
    xp = [xp;xpbl];
    vxp = [vxp;vxpbl];
    wp = [wp;wpbl];
    
    % Update Density
    ni(:) = 0.0;
    ne_prev = ne;
end
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(xp,vxp/3e10,'.','MarkerSize',3);
title('x-Px');
figure;
plot(x,E);
title('Ex');
figure;
xh = (x(1:end-1)+x(2:end))/2;
plot(xh,ne);
title('ne');
figure;
plot(xh,Te);
title('Te');
%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%