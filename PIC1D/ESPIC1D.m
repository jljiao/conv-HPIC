%% Electrostatic Particle-in-Cell (ESPIC) Code (1D)%%%%%%%%%%%%%%%%%%%%%%%%
% The ESPIC1D is an one-dimentional electrostaic PIC code.This code is 
% used to simulate a conterstream plasma system, in which electrostatic 
% collisionless shock can be excited.
% This code is written by J. L. Jiao.
clc;clear;
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In c.g.s Unit
e = 4.8e-10; % Elementary Charge.
me = 9.11e-28; % Electron Mass [g].
mp = 1.67e-24; % Proton Mass [g].
c = 3e10; % Speed of Light [cm/s].
kB = 1.38e-16; % Boltzmann Constant.
%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
t = 5e-12; % Total Simulation Time [s]
dt_factor = 0.8; % Ratio of dt and 1/omega_pi.
Lx = 1e-2; % The Length of the Simulation Box [cm].
NxDebye = 4; % The Cell Numbers in Debye Length.
Npe = 1e6; % The Total Number of the Electron Pesudo-Particles. 
Npi = 1e6; % The Total Number of the Ion Pesudo-Particles. 
Ns = 3; % Smoothing Cell Numbers.

% Plasma Parameters
Z = 1; % Atom Number.
A = 1; % Mass Number.
ne0 = 1.0e20; % Electron Density [/cm^3].
Te0 = 40*1e3; % Electron Temperature [eV].
Ti0 = 1e3; % Ion Temperature [eV].
ux0 = 0.01*c; % Fluid Velocity in X Direction [cm/s].
%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qi = Z*e; % Charge of Ion.
mi = mp*A; % Mass of Ion.
ni0 = ne0/Z; % Density of Ion.
lambda_De = 743*(Te0/ne0)^0.5; % The Debye Length.

dx = lambda_De/NxDebye; % Cell Size.
Nx = floor(Lx/dx); % Cell Numbers.
if(mod(Nx,2)==1)
    Nx = Nx+1;
end
Lx = Nx*dx;
ppce = floor(Npe/Nx); % The Number of the Electron Pesudo-Particles per Cell. 
ppci = floor(Npi/Nx); % The Number of the Ion Pesudo-Particles per Cell. 
dt = dx/c*dt_factor; % Time Step Size.
Nt = floor(t/dt); % Numbers of Time Step.
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids Initialization
x = linspace(-dx,Lx+dx,Nx+3); % Coordinates of the Box.
ni = zeros(1,Nx+3); % Ion Density.
ne = zeros(1,Nx+3); % Electron Density.
E = zeros(1,Nx+3); % Electrostatic Field.
E_average = zeros(1,Nx+3); % Diagnostic Averaged Electrostatic Field.

% Particles Initialization
npe = (Nx/2+1)*ppce;
xpe = [rand(npe,1)*(Lx/2+dx)-dx;...
    rand(npe,1)*(Lx/2+dx)+Lx/2]; % X Position of the Ion Pesudo-Particles.
vxpe = [randn(npe,1)*(2*Te0/5.11e5)^0.5*c+ux0;...
    randn(npe,1)*(2*Te0/5.11e5)^0.5*c-ux0]; % vx of the Ion Pesudo-Particles
wpe = ne0/ppce*ones(2*npe,1); % Weight of the Ion Pesudo-Particles

npi = (Nx/2+1)*ppci;
xpi = [rand(npi,1)*(Lx/2+dx)-dx;...
    rand(npi,1)*(Lx/2+dx)+Lx/2]; % X Position of the Ion Pesudo-Particles.
vxpi = [randn(npi,1)*(2*Ti0/A/9.38e8)^0.5*c+ux0;...
    randn(npi,1)*(2*Ti0/A/9.38e8)^0.5*c-ux0]; % vx of the Ion Pesudo-Particles
wpi = ni0/ppci*ones(2*npi,1); % Weight of the Ion Pesudo-Particles
%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nt
    % Update Density
    ni(:) = 0.0;
    ne(:) = 0.0;
    
    % Calculate Ion Density
    ind0 = floor(xpi/dx);
    ind1 = ind0+1;
    w0 = 1-(xpi/dx-ind0);
    w1 = 1-w0;
    n0 = w0.*wpi;
    n1 = w1.*wpi;
    for j=1:length(xpi)
        ni(ind0(j)+2) = ni(ind0(j)+2)+n0(j);
        ni(ind1(j)+2) = ni(ind1(j)+2)+n1(j);
    end
    ni(1) = ni0;
    ni(end) = ni0;
    ni = movmean(ni,Ns);

    % Calculate Electron Density
    ind0 = floor(xpe/dx);
    ind1 = ind0+1;
    w0 = 1-(xpe/dx-ind0);
    w1 = 1-w0;
    n0 = w0.*wpe;
    n1 = w1.*wpe;
    for j=1:length(xpe)
        ne(ind0(j)+2) = ne(ind0(j)+2)+n0(j);
        ne(ind1(j)+2) = ne(ind1(j)+2)+n1(j);
    end
    ne(1) = ne0;
    ne(end) = ne0;
    ne = movmean(ne,Ns);

    % Calculate Electrostatic Field
    rho = 4*pi*e*(ni-ne);
    phi = poisson1D(rho',0.0,0.0,dx)';
    E(2:end-1) = (phi(3:end)-phi(1:end-2))/dx/2;
    E(1) = 0.0;
    E(end) = 0.0;
%     E = movmean(E,Ns);

    % Push Particles
    vxpi = vxpi+qi*dt/mi*interp1(x,E,xpi);
    xpi = xpi+vxpi*dt;
    vxpe = vxpe-e*dt/me*interp1(x,E,xpe);
    xpe = xpe+vxpe*dt;

    % Ion Injector Boundary for Right Side
    vxpi(xpi>=Lx) = [];
    wpi(xpi>=Lx) = [];
    xpi(xpi>=Lx) = [];
    
    xpbr = Lx+rand(ppci,1)*dx;
    vxpbr = randn(ppci,1)*(2*Ti0/A/9.38e8)^0.5*c-ux0;
    wpbr = ni0/ppci*ones(ppci,1);
    
    xpi = [xpi;xpbr];
    vxpi = [vxpi;vxpbr];
    wpi = [wpi;wpbr];
    
    % Electron Injector Boundary for Right Side
    vxpe(xpe>=Lx) = [];
    wpe(xpe>=Lx) = [];
    xpe(xpe>=Lx) = [];
    
    xpbr = Lx+rand(ppce,1)*dx;
    vxpbr = randn(ppce,1)*(2*Te0/5.11e5)^0.5*c-ux0;
    wpbr = ne0/ppce*ones(ppce,1);
    
    xpe = [xpe;xpbr];
    vxpe = [vxpe;vxpbr];
    wpe = [wpe;wpbr];
    
    % Ion Injector Boundary for Left Side
    vxpi(xpi<=0) = [];
    wpi(xpi<=0) = [];
    xpi(xpi<=0) = [];
    
    xpbl = -rand(ppci,1)*dx;
    vxpbl = randn(ppci,1)*(2*Ti0/A/9.38e8)^0.5*c+ux0;
    wpbl = ni0/ppci*ones(ppci,1);
    
    xpi = [xpi;xpbl];
    vxpi = [vxpi;vxpbl];
    wpi = [wpi;wpbl];
    
    % Electron Injector Boundary for Left Side
    vxpe(xpe<=0) = [];
    wpe(xpe<=0) = [];
    xpe(xpe<=0) = [];
    
    xpbl = -rand(ppce,1)*dx;
    vxpbl = randn(ppce,1)*(2*Te0/5.11e5)^0.5*c+ux0;
    wpbl = ne0/ppce*ones(ppce,1);
    
    xpe = [xpe;xpbl];
    vxpe = [vxpe;vxpbl];
    wpe = [wpe;wpbl];
    
    % Diagnostic Averaged Electrostatic Field
    if(i>Nt-100)
        E_average = E_average+E;
    end
end
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
vxpi = vxpi/3e10;
plot(xpi,vxpi,'.','MarkerSize',3);
title('Ion Phase');

figure;
vxpe = vxpe/3e10;
plot(xpe,vxpe/3e10,'.','MarkerSize',3);
title('Electron Phase');

figure;
E_average = E_average/100;
plot(E_average);
title('Ex');

% figure;
% plot(rho);
% title('\rho');

figure;
plot(ne);
hold on;
plot(ni);
title('Density');
%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%