clear all, close all;
% Parameters
nbrOfPoints = 10000;
r0 = 0.001;
rEnd = 10;
rValues = linspace(r0,rEnd,nbrOfPoints)'; 
h = rValues(2) - rValues(1); % grid spacing

% Wavefunction
f = @(r) 2 * r.* exp(-r);

% Initialize matrix
A = zeros(nbrOfPoints);
for i=2:nbrOfPoints-1 % Inner points
    A(i,i+1) = 1/(h^2);
    A(i,i-1) = 1/(h^2);
    A(i,i) = - 2/(h^2);
end

% Boundary elements
A(1,1) = -2/(h^2);
A(1,2) = 1/(h^2);
A(end,end) = -2/(h^2);
A(end,end-1) = 1/(h^2);

% Hartree potential
Vh = 1./rValues - (1 + 1./rValues).*exp(-2*rValues);

potential1 = -2./rValues + Vh;
potential2 = -1./rValues;
nuclearPotential = potential2;

potF = -0.5*A *f(rValues) + nuclearPotential.*f(rValues);
energyMatrix = potF ./ f(rValues);
mean(energyMatrix(1:end-1))
