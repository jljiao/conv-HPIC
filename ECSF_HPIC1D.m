%% HPIC Code for ECSF 1D Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an one-dimentional (1D) hybrid Particle-in-Cell (HPIC) code. It 
% simulates Electrostatic Charge Separation Field (ECSF) in an 1D steady-
% state hydrogen plasma with step-like density by convolutional method and 
% iterative method. This code is written by Jin-Long Jiao.
% In cgs Unit.
%% Clean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 4.8e-10;                           % Elementary Charge [statcoulomb]
%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = 1;                                 % Atomic Number of Hydrogen
ne_0 = 1e20;                           % Plasma Density [cm^-3]
Te_0 = 40e3;                           % Plasma Temperature [eV]
gamma = 3;                             % Specific Heat Ratio
n_cell_De = 4;                         % Cell Numbers per Debye Length
n_cell = 1000;                         % Total Cell Numbers
err_limit = 1e3;                       % Iterative Absolute Error
n_De_conv = 20;                        % Numbers of Debye Length Contained 
                                       % in Convolutional Kernel Radius
%% Simulation Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_De = 743*(Te_0/ne_0)^0.5;       % Debye Length [cm]
dx = lambda_De/n_cell_De;              % Cell Size [cm]
x = dx*(1:n_cell+1);                   % Grid Coordinates [cm]
ni = [ones(1,n_cell/2), ...
      ones(1,n_cell/2)*0.5]*ne_0;      % Initial Ion Density [cm^-3]
ne = Z*ni;                             % Initial Electron Density [cm^-3]
Te = Te_0*ones(1,n_cell);              % Initial Electron Temperature [eV]
l_De = 743*(Te./ne).^0.5;              % Debye Length on Cell
l_De_c = (l_De(1:end-1)+l_De(2:end))/2;% Debye Length on Grid
n_kernel = n_De_conv*n_cell_De;        % Grid Numbers of Kernel Radius

%% Iterative Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
chi = gamma*l_De_c.^2./(dx^2+ ...
      2*gamma*l_De_c.^2);              % Coefficient \chi
E_iter = zeros(1,n_cell+1);            % Electrostatic Field Ex

% Iteration
err = err_limit;
while(err >= err_limit)
    E0 = E_iter;
    E_iter(2:end-1) = chi.*((E0(1:end-2)+E0(3:end))- ...
                            4*pi*e*dx*(ne(2:end)-ne(1:end-1)));
    E_iter(1) = 0.0; E_iter(end) = 0.0;
    err = max(abs(E_iter-E0));
end

%% Convolutional Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
l_De_conv = [l_De(1)*ones(1,n_kernel),l_De,l_De(end)*ones(1,n_kernel)];
E_conv = zeros(1,n_cell+1);            % Electrostatic Field Ex
ne_conv = [ne(1)*ones(1,n_kernel),ne,ne(end)*ones(1,n_kernel)];
kernel = zeros(1,1+2*n_kernel);        % Convolutional Kernel

% Convlution
for j=1:n_cell
    kernel(n_kernel+1) = 1;
    for i=1:n_kernel
        kernel(n_kernel+1+i) = kernel(n_kernel+i)* ...
                               exp(-dx/l_De_conv(j+n_kernel+i)/gamma^0.5);
        kernel(n_kernel+1-i) = kernel(n_kernel+2-i)* ...
                               exp(-dx/l_De_conv(j+n_kernel-i)/gamma^0.5);
    end
    kernel = kernel./l_De_conv(j:j+2*n_kernel);
    kernel = kernel/sum(kernel);
    ne(j) = sum(ne_conv(j:j+2*n_kernel).*kernel);
end

E_conv(2:end-1) = -4*pi*e*gamma*((l_De(2:end)+l_De(1:end-1))/2).^2.* ...
                  (ne(2:end)-ne(1:end-1))/dx;
E_conv(1) = 0.0; E_conv(end) = 0.0;    % Open Boundary for Ex
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Figure
figure;
x0 = [0,0,500*dx*1e4,500*dx*1e4,1000*dx*1e4,1000*dx*1e4];
y0 = [0,1,1,0.5,0.5,0];
fill(x0,y0,'g');
hold on;
plot(x*1e4,E_iter*3e4/1e11,'-b','LineWidth',2);
hold on;
plot(x*1e4,E_conv*3e4/1e11,'--r','LineWidth',2);

% Plot Setting
xlim([0,n_cell*dx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$E_x~[10^{11}\rm{V/m}]$', 'interpreter', 'latex');
legend('Density [$10^{21}$cm$^{-3}$]','Iterative Method', ...
       'Convolutional Method','interpreter','latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);
%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%