%% PIC Code for CES 1D1V Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an one-dimentional for space and velocity (1D1V) fully kinetic 
% Particle-in-Cell (PIC) code. It simulates Collisionless Electrostatic 
% Shock (CES) in an 1D supersonic counterstreaming plasmas system. 
% This code is written by Jin-Long Jiao.
% In cgs Unit.
%% Clean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
tic;                                   % Measure Elapsed Time Start
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 4.8e-10;                           % Elementary Charge [statcoulomb]
m_e = 9.11e-28;                        % Electron Mass [g]
m_p = 1.67e-24;                        % Proton Mass [g]
c = 3e10;                              % Speed of Light [cm/s]

%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
t = 5e-12;                             % Total Simulation Time [s]
dt_factor = 0.1;                       % Factor of Time Step
Lx = 1e-2;                             % Simulation Box Length [cm]
n_cell_De = 4;                         % Cell Numbers per Debye Length
n_part_e = 1e6;                        % Numbers of Electron MacroParticles
n_part_i = 1e6;                        % Numbers of Ion Macro-Particles

% Plasma Parameters
Z = 1;                                 % Ion Atom Number
A = 1;                                 % Ion Mass Number
ne_0 = 1.0e20;                         % Electron Density [cm^-3]
Te_0 = 40*1e3;                         % Electron Temperature [eV]
Ti_0 = 1e3;                            % Ion Temperature [eV]
ux_0 = 0.01*c;                         % Plasma Colliding Velocity [cm/s]

%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qi = Z*e;                              % Ion Charge
mi = m_p*A;                            % Ion Mass
ni_0 = ne_0/Z;                         % Ion Initial Density [cm^-3]
omega_pe = 5.64e4*(Z*ni_0)^0.5;        % Electron Plasma Frequency [s^-1]
lambda_De = 743*(Te_0/ne_0)^0.5;       % Debye Length [cm]

dx = lambda_De/n_cell_De;              % Cell Size [cm]
n_cell = floor(Lx/dx);                 % Total Cell Numbers
if(mod(n_cell,2)==1)
    n_cell = n_cell+1;                 % Reset Total Cell Numbers
end
Lx = n_cell*dx;                        % Reset Simulation Box Length [cm]
ppc_e = floor(n_part_e/n_cell);        % Electron Macro-Particles per Cell 
ppc_i = floor(n_part_i/n_cell);        % Ion Macro-Particles per Cell 
dt = 1/omega_pe*dt_factor;             % Time Step Size [s]
n_t = floor(t/dt);                     % Total Time Step Numbers

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids Initialization
x = linspace(-dx,Lx+dx,n_cell+3);      % Grid Coordinates [cm]
ni = zeros(1,n_cell+3);                % Define Ion Density
ne = zeros(1,n_cell+3);                % Define Electron Density
Ex = zeros(1,n_cell+3);                % Define Electrostatic Field
Ex_aver = zeros(1,n_cell+3);           % Define Averaged Ex

% Particles Initialization
n_p_e = (n_cell/2+1)*ppc_e;            % Electron Numbers
x_p_e = [rand(n_p_e,1)*(Lx/2+dx)-dx;...
    rand(n_p_e,1)*(Lx/2+dx)+Lx/2];     % Electron Coordinates [cm]
vx_p_e = [randn(n_p_e,1)*(2*Te_0/5.11e5)^0.5*c+ux_0;randn(n_p_e,1)* ...
          (2*Te_0/5.11e5)^0.5*c-ux_0]; % Electron Velocities [cm/s]
w_e = ne_0/ppc_e*ones(2*n_p_e,1);      % Electron Weights

n_p_i = (n_cell/2+1)*ppc_i;            % Ion Numbers
x_p_i = [rand(n_p_i,1)*(Lx/2+dx)-dx;...
    rand(n_p_i,1)*(Lx/2+dx)+Lx/2];     % Ion Coordinates [cm]
vx_p_i = [randn(n_p_i,1)*(2*Ti_0/A/9.38e8)^0.5*c+ux_0;randn(n_p_i,1)* ...
         (2*Ti_0/A/9.38e8)^0.5*c-ux_0];% Ion Velocities [cm/s]
w_i = ni_0/ppc_i*ones(2*n_p_i,1);      % Ion Weights

%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n_t
    % Density Initialization
    ni(:) = 0.0;
    ne(:) = 0.0;

    % Calculate Ion Density
    ind = floor(x_p_i/dx);
    w0 = 1-(x_p_i/dx-ind);
    w1 = 1-w0;
    n0 = w0.*w_i;
    n1 = w1.*w_i;
    ind = ind+2;
    for j=1:length(x_p_i)
        ni(ind(j)) = ni(ind(j))+n0(j);
        ni(ind(j)+1) = ni(ind(j)+1)+n1(j);
    end
    ni(1) = ni_0; ni(end) = ni_0;      % Injector Boundary for Density
    ni = movmean(ni,3);                % Density Smoothing

    % Calculate Electron Density
    ind = floor(x_p_e/dx);
    w0 = 1-(x_p_e/dx-ind);
    w1 = 1-w0;
    n0 = w0.*w_e;
    n1 = w1.*w_e;
    ind = ind+2;
    for j = 1:length(x_p_e)
        ne(ind(j)) = ne(ind(j))+n0(j);
        ne(ind(j)+1) = ne(ind(j)+1)+n1(j);
    end
    ne(1) = ne_0; ne(end) = ne_0;      % Injector Boundary for Density
    ne = movmean(ne,3);                % Density Smoothing

    % Calculate Electrostatic Field
    rho = 4*pi*e*(ne-ni);
    phi = poisson1D(rho',0.0,0.0,dx)';
    Ex(2:end-1) = -(phi(3:end)-phi(1:end-2))/dx/2;
    Ex(1) = 0.0; Ex(end) = 0.0;       % Open Boundary for Field

    % Push Particles
    vx_p_i = vx_p_i+qi*dt/mi*interp1(x,Ex,x_p_i);
    x_p_i = x_p_i+vx_p_i*dt;
    vx_p_e = vx_p_e-e*dt/m_e*interp1(x,Ex,x_p_e);
    x_p_e = x_p_e+vx_p_e*dt;

    % Ion Injector Boundary for Right Side
    vx_p_i(x_p_i>=Lx) = [];
    w_i(x_p_i>=Lx) = [];
    x_p_i(x_p_i>=Lx) = [];
    
    x_p_br = Lx+rand(ppc_i,1)*dx;
    vx_p_br = randn(ppc_i,1)*(2*Ti_0/A/9.38e8)^0.5*c-ux_0;
    w_br = ni_0/ppc_i*ones(ppc_i,1);
    
    x_p_i = [x_p_i;x_p_br];
    vx_p_i = [vx_p_i;vx_p_br];
    w_i = [w_i;w_br];
    
    % Electron Injector Boundary for Right Side
    vx_p_e(x_p_e>=Lx) = [];
    w_e(x_p_e>=Lx) = [];
    x_p_e(x_p_e>=Lx) = [];
    
    x_p_br = Lx+rand(ppc_e,1)*dx;
    vx_p_br = randn(ppc_e,1)*(2*Te_0/5.11e5)^0.5*c-ux_0;
    w_br = ne_0/ppc_e*ones(ppc_e,1);
    
    x_p_e = [x_p_e;x_p_br];
    vx_p_e = [vx_p_e;vx_p_br];
    w_e = [w_e;w_br];
    
    % Ion Injector Boundary for Left Side
    vx_p_i(x_p_i<=0) = [];
    w_i(x_p_i<=0) = [];
    x_p_i(x_p_i<=0) = [];
    
    x_p_bl = -rand(ppc_i,1)*dx;
    vx_p_bl = randn(ppc_i,1)*(2*Ti_0/A/9.38e8)^0.5*c+ux_0;
    w_bl = ni_0/ppc_i*ones(ppc_i,1);
    
    x_p_i = [x_p_i;x_p_bl];
    vx_p_i = [vx_p_i;vx_p_bl];
    w_i = [w_i;w_bl];
    
    % Electron Injector Boundary for Left Side
    vx_p_e(x_p_e<=0) = [];
    w_e(x_p_e<=0) = [];
    x_p_e(x_p_e<=0) = [];
    
    x_p_bl = -rand(ppc_e,1)*dx;
    vx_p_bl = randn(ppc_e,1)*(2*Te_0/5.11e5)^0.5*c+ux_0;
    w_bl = ne_0/ppc_e*ones(ppc_e,1);
    
    x_p_e = [x_p_e;x_p_bl];
    vx_p_e = [vx_p_e;vx_p_bl];
    w_e = [w_e;w_bl];
    
    % Diagnostic Averaged Electrostatic Field
    if(i > n_t-100)
        Ex_aver = Ex_aver+Ex;
    end
end
toc;
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(x_p_i*1e4,vx_p_i/3e10,'.','MarkerSize',3);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('Ion $v_x~[c]$', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

% figure;
% plot(x_p_e*1e4,vx_p_e/3e10,'.','MarkerSize',3);
% xlim([0,Lx*1e4]);
% xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
% ylabel('Electron $v_x~[c]$', 'interpreter', 'latex');
% 
% set(gcf,'Position',[100 100 800 500]);
% set(gca,'Position',[.15 .2 .7 .7]);
% set(gcf,'color',[1,1,1]);
% set(gca,'Fontsize',18);

figure;
plot(x*1e4,Ex_aver*3e4/1e12);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$E_x$ [TV/m]', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

figure;
plot(x*1e4,ne/1e20);
hold on;
plot(x*1e4,ni/1e20);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('Density [10$^{20}$cm$^{-3}$]', 'interpreter', 'latex');
legend('$n_e$','$n_i$', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

%% Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB function is written by Koorosh Gobal. Jin-Long Jiao gets this
% function from MATLAB website: https://www.math-works.com/matlabcentral/
% fileexchange/49801-poisson1d-f-uleft-uright-n-l

% POISSON1D solves the poisson equation d2U/dX2 = F.
% u = poisson1D(f,Uleft,Uright,N,L)
% f: a vector representing the right-hand-side
% Uleft: dirichlet boundary condition at u(0)
% Uright: dirichlet boundary condition at u(L)
% N: Number of nodes
% L: Length of domain
function out = poisson1D(f,Uleft,Uright,dx)
    N = length(f);
    uB = zeros(length(f),1);
    uB(2) = Uleft;
    uB(end-1) = Uright;
    f = f*dx^2-uB;
    b = dst(f);
    m = (1:length(b))';
    a = b./(2*(cos(m*pi/(N-1))-1));
    uSOL = idst(a);
    uSOL(1) = Uleft;
    uSOL(end) = Uright;
    out = uSOL;
end
%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%