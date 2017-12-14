clc, clear all, close all;
% Parameters
nbrOfPoints = 10;
rEnd = 10;
h = rEnd/nbrOfPoints; % grid spacing
r0 = h;
rValues = linspace(r0,rEnd,nbrOfPoints)';

% Initial guess for wavefunction
f = exp(-rValues)/sqrt(pi);

% Initialize matrix
A = zeros(nbrOfPoints);
for i=2:nbrOfPoints-1 % Inner points
    A(i,i+1) = 1/(h^2);
    A(i,i-1) = 1/(h^2);
    A(i,i) = - 2/(h^2);
end

% Edges
A(1,1) = -2/(h^2);
A(1,2) = 1/(h^2);
A(end,end) = -2/(h^2);
A(end,end-1) = 1/(h^2);

% Boundary conditions
c = zeros(nbrOfPoints,1);
c(end) = - 1/(h^2);

energy = [0,10000];
threshold = 0.001;
%while(abs(energy(end)-energy(end-1))>threshold)
   % Compute electron density
   ns = f.^2;
    % Compute Hartree potential
    ns2 = ns.*rValues;
    U = A\(- 4*pi*ns2 + c);
    Vh = U./rValues;
    
    potential = -2 * diag(1./rValues) + diag(Vh);
    H = -0.5*A + potential;
    
    [eigVec,eigVal] = eig(H);
    epsilon = L(1,1)
%end
