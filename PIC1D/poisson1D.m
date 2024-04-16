% This MATLAB function is written by Koorosh Gobal.
% J.L. Jiao gets this file from the website of MATLAB:(https://www.math-
% works.com/matlabcentral/fileexchange/49801-poisson1d-f-uleft-uright-n-l)

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