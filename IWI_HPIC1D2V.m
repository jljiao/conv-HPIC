%% HPIC Code for IWI 1D2V Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an one-dimentional for space and two-dimentional for velocity 
% (1D2V) hybrid Particle-in-Cell (HPIC) code. It simulates Ion-Weibel 
% Instability (IWI) in an 1D supersonic counterstreaming plasmas system by 
% convolutional method and iterative method. 
% This code is written by Jin-Long Jiao.
% In cgs Unit.
%% Clean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
%% Fundamental Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 4.8e-10;                           % Elementary Charge [statcoulomb]
m_p = 1.67e-24;                        % Proton Mass [g]
m_e = 9.11e-28;                        % Electron Mass [g]
c = 3e10;                              % Speed of Light [cm/s]
%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Parameters
t = 20e-12;                            % Total Simulation Time [s]
dt_factor = 0.1;                       % Factor of Time Step
Lx = 2e-3;                             % Simulation Box Length [cm]
n_cell_De = 4;                         % Cell Numbers per Debye Length
n_part = 2e6;                          % Numbers of Ion Macro-Particles
eta = 2;                               % Correction Factor for ICS Eq.

% Plasma Parameters
Z = 1;                                 % Ion Atom Number
A = 1;                                 % Ion Mass Number
cmr = Z*e/(A*m_p);                     % Charge-to-Mass Ratio of Ion
ne_0 = 2.0e20;                         % Electron Density [cm^-3]
Te_0 = 3.5e3;                          % Electron Temperature [eV]
Ti_0 = 1e2;                            % Ion Temperature [eV]
uy_0 = 0.1*c;                          % Plasma Colliding Velocity [cm/s]
gamma = 2;                             % Specific Heat Ratio

% Physics Option
phys_opt_Bz = 2;                       % 1: Iterative Method for Bz
                                       % 2: Convolutional Method for Bz
switch  phys_opt_Bz
    case 1
        err_Bz_limit = 1e0;            % Iterative Absolute Error for Bz
    case 2
        n_De_conv_Bz = 100;            % Numbers of Debye Length Contained 
                                       % in Convolutional Kernel Radius
        n_kernel_Bz = n_De_conv_Bz* ...
                      n_cell_De;       % Grid Numbers of Kernel Radius
        kernel_Bz = ones(1,1+2* ...
                         n_kernel_Bz); % Convolutional Kernel Function
        alpha = eta*m_e*c/e;
end

phys_opt_Ex = 2;                       % 1: Iterative Method for Ex.
                                       % 2: Convolutional Method for Ex.
switch  phys_opt_Ex
    case 1
        err_Ex_limit = 1e0;            % Iterative Absolute Error for Ex
    case 2
        n_De_conv_Ex = 20;             % Numbers of Debye Length Contained 
                                       % in Convolutional Kernel Radius
        n_kernel_ne = n_De_conv_Ex* ...
                      n_cell_De;       % Grid Numbers of Kernel Radius
        kernel_ne = ones(1,1+2* ...
                         n_kernel_ne); % Convolutional Kernel Function
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
n_t_diag = 10;                         % Diagnostic Time Step Numbers
Bz_diag = zeros(floor(n_t/n_t_diag), ...
                n_cell+4);             % Diagnostic Magnetic Field [G]
ind_diag = 1;                          % Diagnostic Index

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids Initialization
x = linspace(-1.5*dx,Lx+1.5*dx, ...
             n_cell+4);                % Grid Coordinates [cm]
ni = zeros(1,n_cell+3);                % Define Ion Density
uy = zeros(1,n_cell+3);                % Define Ion Fluid Velocity
Ex = zeros(1,n_cell+4);                % Define Electrostatic Field
Bz = zeros(1,n_cell+4);                % Define Magnetic Field
F = zeros(1,n_cell+4);                 % Define Temporary variable

% Particles Initialization
n_p = (n_cell+2)*ppc;                  % Particle Numbers of Each Plasma
x_p = rand(2*n_p,1)*(Lx+2*dx)-dx;      % Particle Coordinates [cm]
vx_p = [randn(n_p,1)*v_th; ...
        randn(n_p,1)*v_th];            % Particle Velocities in X [cm/s]
vy_p = [randn(n_p,1)*v_th+uy_0; ...
        randn(n_p,1)*v_th-uy_0];       % Particle Velocities in Y [cm/s]
w = 0.5*ni_0/ppc*ones(2*n_p,1);        % Particle Weights

%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n_t
    % Density Initialization
    ni(:) = 0.0;
    % Fluid Velocity Initialization
    uy(:) = 0.0;

    % Calculate Density and Fluid Velocity
    ind = floor(x_p/dx);
    w0 = 1-(x_p/dx-ind);
    w1 = 1-w0;
    n0 = w0.*w;
    n1 = w1.*w;
    u0 = w0.*w.*vy_p;
    u1 = w1.*w.*vy_p;
    ind = ind+2;
    for j=1:length(x_p)
        ni(ind(j)) = ni(ind(j))+n0(j);
        ni(ind(j)+1) = ni(ind(j)+1)+n1(j);
        uy(ind(j)) = uy(ind(j))+u0(j);
        uy(ind(j)+1) = uy(ind(j)+1)+u1(j);
    end
    ni(1) = ni(end-1); ni(end) = ni(2);% Periodic Boundary for Density
    uy(1) = uy(end-1); uy(end) = uy(2);% Periodic Boundary for Velocity
    uy = uy./ni;

    ni = movmean(ni,3);                % Density Smoothing
    uy = movmean(uy,3);                % Velocity Smoothing
    ne = Z*ni;                         % Estimated Electron Density

    % Calculate Magnetic Field (Bz)
    ls = ((m_e*c^2/4/pi/e^2)./ne).^0.5;% Plasma Skin Depth on Cell
    ls_c = (ls(1:end-1)+ls(2:end))/2;  % Plasma Skin Depth on Grid
    switch  phys_opt_Bz
        case 1
            ls2 = ls_c.^2;             % Square of Plasma Skin Depth
            delta = ls2./(2*ls2+dx^2/eta);
            ksi = m_e*c/e./ls2*dx;
 
            err = err_Bz_limit;
            while(err >= err_Bz_limit)
                B0 = Bz;
                Bz(2:end-1) = delta.*(B0(1:end-2)+B0(3:end)+ ...
                                      ksi.*(uy(2:end)-uy(1:end-1)));
                Bz(1) = Bz(end-1);
                Bz(end) = Bz(2);       % Periodic Boundary for Bz
                err = max(abs(Bz-B0));
            end
        case 2
            ls_conv = [ls_c(end-n_kernel_Bz-1:end-1),ls_c, ...
                       ls_c(2:n_kernel_Bz+2)];
            curl_u = (uy(2:end)-uy(1:end-1))/dx;
            curl_u_conv = [curl_u(end-n_kernel_Bz-1:end-1),curl_u, ...
                           curl_u(2:n_kernel_Bz+2)];

            for j = 1:n_cell+4
                kernel_Bz(n_kernel_Bz+1) = 1;
                for k = 1:n_kernel_Bz
                    kernel_Bz(n_kernel_Bz+1+k) = ...
                                              kernel_Bz(n_kernel_Bz+k)* ...
                                 exp(-dx/ls_conv(j+n_kernel_Bz+k)/eta^0.5);
                    kernel_Bz(n_kernel_Bz+1-k) = ...
                                            kernel_Bz(n_kernel_Bz+2-k)* ...
                                 exp(-dx/ls_conv(j+n_kernel_Bz-k)/eta^0.5);
                end
                kernel_Bz = kernel_Bz/sum(kernel_Bz);
                Bz(j) = alpha*sum(curl_u_conv(j:j+2*n_kernel_Bz).* ...
                                  kernel_Bz);
            end
    end

    % Calculate Electrostatic Field
    Te = Te_0*(ne/ne_0).^(gamma-1);    % EOS with Adiabatic Approximation
    l_De = 743*(Te./ne).^0.5;          % Debye Length
    switch  phys_opt_Ex
        case 1
            l_De_c = (l_De(1:end-1)+ ...
                      l_De(2:end))/2;  % Debye Length on Grid
            chi = gamma*l_De_c.^2./(dx^2+2*gamma*l_De_c.^2);
            kappa = dx^2/gamma./l_De_c.^2;

            err = err_Ex_limit;
            while(err >= err_Ex_limit)
                E0 = Ex;
                F0 = (uy(1:end-1)+uy(2:end)).*Bz(2:end-1)/2/c+ ...
                     (Bz(3:end).^2-Bz(1:end-2).^2)./(ne(1:end-1)+ ...
                      ne(2:end))/8/pi/e/dx;
                Ex(2:end-1) = chi.*(E0(1:end-2)+E0(3:end)-4*pi*e*dx*Z* ...
                              (ni(2:end)-ni(1:end-1))-kappa.*F0);
                Ex(1) = Ex(end-1);
                Ex(end) = Ex(2);       % Periodic Boundary for Ex
                err = max(abs(Ex-E0));
            end
        case 2
            ne_conv = [ne(end-1-n_kernel_ne:end-2),ne,ne(3:n_kernel_ne+2)];
            l_De_conv = [l_De(end-1-n_kernel_ne:end-2),l_De, ...
                        l_De(3:n_kernel_ne+2)];

            F(2:end-1) = (uy(1:end-1)+uy(2:end)).*Bz(2:end-1)/2/c+ ...
                         (Bz(3:end).^2-Bz(1:end-2).^2)./ ...
                         (ne(1:end-1)+ne(2:end))/8/pi/e/dx;
            F(1) = F(end-1); F(end) = F(2);
            div_F = (F(2:end)-F(1:end-1))/dx/(4*pi*e);
            div_F_conv = [div_F(end-1-n_kernel_ne:end-2),div_F, ...
                          div_F(3:n_kernel_ne+2)];

            for j = 1:n_cell+3
                kernel_ne(n_kernel_ne+1) = 1;
                for k = 1:n_kernel_ne
                    kernel_ne(n_kernel_ne+1+k) = ...
                                              kernel_ne(n_kernel_ne+k)* ...
                              exp(-dx/l_De_conv(j+n_kernel_ne+k)/gamma^0.5);
                    kernel_ne(n_kernel_ne+1-k) = ...
                                           kernel_ne(n_kernel_ne+2-k)* ...
                              exp(-dx/l_De_conv(j+n_kernel_ne-k)/gamma^0.5);
                end

                kernel_ne = kernel_ne./l_De_conv(j:j+2*n_kernel_ne);
                kernel_ne = kernel_ne/sum(kernel_ne);
                ne(j) = sum((ne_conv(j:j+2*n_kernel_ne)+ ...
                            div_F_conv(j:j+2*n_kernel_ne)).*kernel_ne);
            end

            Te = Te_0*(ne/ne_0).^ ...
                 (gamma-1);            % EOS with Adiabatic Approximation
            l_De = 743*(Te./ne).^0.5;   % Debye Length
            F(2:end-1) = (uy(1:end-1)+uy(2:end)).*Bz(2:end-1)/2/c+ ...
                         (Bz(3:end).^2-Bz(1:end-2).^2)./ ...
                         (ne(1:end-1)+ne(2:end))/8/pi/e/dx;
            Ex(2:end-1) = -4*pi*e*gamma*((l_De(2:end)+ ...
                          l_De(1:end-1))/2).^2.*...
                          (ne(2:end)-ne(1:end-1))/dx-F(2:end-1);
            Ex(1) = Ex(end-1); 
            Ex(end) = Ex(2);           % Periodic Boundary for Ex
    end

    % Push Ions
    Bz_p = interp1(x,Bz,x_p);
    Ex_p = interp1(x,Ex,x_p);
    vx_p = vx_p + cmr*dt*(Ex_p+vy_p/c.*Bz_p);
    vy_p = vy_p - cmr*dt*vx_p/c.*Bz_p;
    x_p = x_p + dt*vx_p;

    % Periodic Boundary
    x_p(x_p>Lx+dx) = x_p(x_p>Lx+dx) - (Lx+2*dx);
    x_p(x_p<-dx) = x_p(x_p<-dx) + (Lx+2*dx);

    % Magnetic Field Diagnostic
    if(mod(i,n_t_diag)==0)
        Bz_diag(ind_diag,:) = Bz;
        ind_diag = ind_diag+1;
    end
end
%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Figure
figure;
plot(x_p*1e4,vy_p/3e10,'.','MarkerSize',3);
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$v_x~[c]$', 'interpreter', 'latex');

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

figure;
Bm = max(max(abs(Bz_diag*1e-6)));
imagesc(x*1e4,linspace(0,t,ind_diag)*1e12,Bz_diag*1e-6);
axis xy;
xlim([0,Lx*1e4]);
xlabel('$x~[\mu\rm{m}]$', 'interpreter', 'latex');
ylabel('$t$ [ps]', 'interpreter', 'latex');
cb = colorbar;
clim([-Bm,Bm]);
cb.Title.String = 'B_z [100T]';

set(gcf,'Position',[100 100 800 500]);
set(gca,'Position',[.15 .2 .7 .7]);
set(gcf,'color',[1,1,1]);
set(gca,'Fontsize',18);

%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%