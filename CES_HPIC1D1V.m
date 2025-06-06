%% HPIC Code for CES 1D1V Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an one-dimentional for space and velocity (1D1V) hybrid Particle-
% in-Cell (HPIC) code. It simulates Collisionless Electrostatic Shock (CES)
% in an 1D supersonic counterstreaming plasmas system by convolutional 
% method and iterative method. This code is written by Jin-Long Jiao.
% In cgs Unit.
%% Clean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
tic;                                   % Measure Elapsed Time Start
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 4.8e-10;                           % Elementary Charge [statcoulomb]
m_p = 1.67e-24;                        % Proton Mass [g]
c = 3e10;                              % Speed of Light [cm/s]
%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
t = 5e-12;                             % Total Simulation Time [s]
dt_factor = 0.1;                       % Factor of Time Step
Lx = 1e-2;                             % Simulation Box Length [cm]
n_cell_De = 4;                         % Cell Numbers per Debye Length
n_part = 1e6;                          % Numbers of Ion Macro-Particles 

% Plasma Parameters
Z = 1;                                 % Ion Atom Number
A = 1;                                 % Ion Mass Number
ne_0 = 1.0e20;                         % Electron Density [cm^-3]
Te_0 = 40e3;                           % Electron Temperature [eV]
Ti_0 = 1e3;                            % Ion Temperature [eV]
ux_0 = 0.01*c;                         % Plasma Colliding Velocity [cm/s]
gamma = 3;                             % Specific Heat Ratio

% Physics Option
phys_opt = 1;                          % 1: Iterative Method
                                       % 2: Convolutional Method
switch  phys_opt
    case 1
        err_limit = 1e3;               % Iterative Absolute Error
    case 2
        n_De_conv = 10;                % Numbers of Debye Length Contained 
                                       % in Convolutional Kernel Radius
        n_kernel = n_De_conv*n_cell_De;% Grid Numbers of Kernel Radius
        kernel = ones(1,1+2*n_kernel); % Convolutional Kernel Function
end
%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qi = Z*e;                              % Ion Charge
mi = m_p*A;                            % Ion Mass
ni_0 = ne_0/Z;                         % Ion Initial Density [cm^-3]
lambda_De = 743*(Te_0/ne_0)^0.5;       % Debye Length [cm]
omega_pi = 1.32e3*Z*(ni_0/A)^0.5;      % Ion Plasma Frequency [s^-1]
v_th = (2*Ti_0/A/9.38e8)^0.5*c;        % Ion Thermal Velocity [cm/s]

dx = lambda_De/n_cell_De;              % Cell Size [cm]
n_cell = floor(Lx/dx);                 % Total Cell Numbers
if(mod(n_cell,2)==1)
    n_cell = n_cell+1;                 % Reset Total Cell Numbers
end
Lx = n_cell*dx;                        % Reset Simulation Box Length [cm]
ppc = floor(n_part/n_cell);            % Ion Macro-Particles per Cell
dt = 1/omega_pi*dt_factor;             % Time Step Size [s]
n_t = floor(t/dt);                     % Total Time Step Numbers
%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids Initialization
x = linspace(-1.5*dx,Lx+1.5*dx, ...
             n_cell+4);                % Grid Coordinates [cm]
ni = zeros(1,n_cell+3);                % Define Ion Density
Ex = zeros(1,n_cell+4);                % Define Electrostatic Field

% Particles Initialization
n_p = (n_cell/2+1)*ppc;                % Particle Numbers of Each Plasma
x_p = [rand(n_p,1)*(Lx/2+dx)-dx; ...
       rand(n_p,1)*(Lx/2+dx)+Lx/2];    % Particle Coordinates [cm]
vx_p = [randn(n_p,1)*v_th+ux_0; ...
        randn(n_p,1)*v_th-ux_0];       % Particle Velocities [cm/s]
w = ni_0/ppc*ones(2*n_p,1);            % Particle Weights
%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n_t
    % Density Initialization
    ni(:) = 0.0;

    % Calculate Densities
    ind = floor(x_p/dx);
    w0 = 1-(x_p/dx-ind);
    w1 = 1-w0;
    n0 = w0.*w;
    n1 = w1.*w;
    ind = ind+2;
    for j=1:length(x_p)
        ni(ind(j)) = ni(ind(j))+n0(j);
        ni(ind(j)+1) = ni(ind(j)+1)+n1(j);
    end
    ni(1) = ni_0; ni(end) = ni_0;      % Injector Boundary for Density
    ni = movmean(ni,3);                % Density Smoothing
    ne = Z*ni;                         % Estimated Electron Density
    
    % Calculate Electron Temperature
    Te = Te_0*(ne/ne_0).^(gamma-1);    % EOS with Adiabatic Approximation
    % Calculate Debye Length
    l_De = 743*(Te./ne).^0.5;          % Debye Length on Cell

    % Calculate Electrostatic Field (Ex)
    switch  phys_opt
        case 1                         % Iterative Method
            l_De_c = (l_De(1:end-1)+ ...
                      l_De(2:end))/2;  % Debye Length on Grid
            chi = gamma*l_De_c.^2./(dx^2+ ...
                  2*gamma*l_De_c.^2);  % Coefficient \chi

            err = err_limit;
            while(err>=err_limit)
                E0 = Ex;
                Ex(2:end-1) = chi.*((E0(1:end-2)+E0(3:end))- ...
                                    4*pi*e*dx*Z*(ni(2:end)-ni(1:end-1)));
                Ex(1) = 0.0;
                Ex(end) = 0.0;         % Open Boundary for Ex
                err = max(abs(Ex-E0));
            end
        case 2                         % Convolutional Method
            ne_conv = [ne_0*ones(1,n_kernel),ne,ne_0*ones(1,n_kernel)];
            l_De_conv = [l_De(1)*ones(1,n_kernel),l_De, ...
                         l_De(end)*ones(1,n_kernel)];
            for j = 1:n_cell+3
                kernel(n_kernel+1) = 1;
                for k = 1:n_kernel
                    kernel(n_kernel+1+k) = kernel(n_kernel+k)* ...
                                exp(-dx/l_De_conv(j+n_kernel+k)/gamma^0.5);
                    kernel(n_kernel+1-k) = kernel(n_kernel+2-k)* ...
                                exp(-dx/l_De_conv(j+n_kernel-k)/gamma^0.5);
                end

                kernel = kernel./l_De_conv(j:j+2*n_kernel);
                kernel = kernel/sum(kernel);
                ne(j) = sum(ne_conv(j:j+2*n_kernel).*kernel);
            end

            Te = Te_0*(ne/ne_0).^ ...
                 (gamma-1);            % EOS with Adiabatic Approximation
            l_De = 743*(Te./ne).^0.5;  % Debye Length on Cell
            l_De_c = (l_De(1:end-1)+ ...
                      l_De(2:end))/2;  % Debye Length on Grid
            Ex(2:end-1) = -4*pi*e*gamma*l_De_c.^2.*...
                          (ne(2:end)-ne(1:end-1))/dx;
            Ex(1) = 0.0;
            Ex(end) = 0.0;             % Open Boundary for Ex
    end

    % Push Ions
    vx_p = vx_p+qi*dt/mi*interp1(x,Ex,x_p);
    x_p = x_p+vx_p*dt;
    
    % Injector Boundary for x+
    vx_p(x_p>=Lx) = [];
    w(x_p>=Lx) = [];
    x_p(x_p>=Lx) = [];
    
    x_p_br = Lx+rand(ppc,1)*dx;
    vx_p_br = randn(ppc,1)*v_th-ux_0;
    w_br = ni_0/ppc*ones(ppc,1);
    
    x_p = [x_p;x_p_br];
    vx_p = [vx_p;vx_p_br];
    w = [w;w_br];
    
    % Injector Boundary for x-
    vx_p(x_p<=0) = [];
    w(x_p<=0) = [];
    x_p(x_p<=0) = [];
    
    x_p_bl = -rand(ppc,1)*dx;
    vx_p_bl = randn(ppc,1)*v_th+ux_0;
    w_bl = ni_0/ppc*ones(ppc,1);
    
    x_p = [x_p;x_p_bl];
    vx_p = [vx_p;vx_p_bl];
    w = [w;w_bl];
end

toc;                                   % Measure Elapsed Time End
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Figure
figure;
plot(x_p*1e4,vx_p/3e10,'.','MarkerSize',3);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$v_x~[c]$', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

figure;
plot(x*1e4,Ex*3e4/1e12);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$E_x$ [TV/m]', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

% xc = (x(1:end-1)+x(2:end))/2;          % Cell Center Coordinates
% figure;
% plot(xc*1e4,ne);
% xlim([0,Lx*1e4]);
% xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
% ylabel('$n_e$ [cm$^{-3}$]', 'interpreter', 'latex');
% 
% set(gcf,'Position',[100 100 800 500]);
% set(gca,'Position',[.15 .2 .7 .7]);
% set(gcf,'color',[1,1,1]);
% set(gca,'Fontsize',18);
% 
% figure;
% plot(xc*1e4,Te);
% xlim([0,Lx*1e4]);
% xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
% ylabel('$T_e$ [eV]', 'interpreter', 'latex');
% 
% set(gcf,'Position',[100 100 800 500]);
% set(gca,'Position',[.15 .2 .7 .7]);
% set(gcf,'color',[1,1,1]);
% set(gca,'Fontsize',18);
%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%