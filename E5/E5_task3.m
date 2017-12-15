clc, clear all, close all;
% Parameters
nbrOfPoints = 1000;
rEnd = 10;
h = rEnd/nbrOfPoints; % grid spacing
r0 = h;
rValues = linspace(r0,rEnd,nbrOfPoints)';

% Make initial guess for wavefunction
f_init = @(r) exp(-r)/sqrt(pi);

f = f_init(rValues);

% Compute electron density
ns = f.^2;

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
threshold = 0.000001;
i = 1;
% Self-consistency loop
while(abs(energy(end)-energy(end-1))>threshold)
    ns = f.^2;
    % Compute Hartree potential
    ns2 = ns.*rValues;
    U = A\(- 4*pi*ns2 + c);
    Vh = U./rValues;
    % Compute Hamiltonian
    potential = -2 * diag(1./rValues) + diag(Vh);
    H = -0.5*A + potential;
    [K,L] = eig(H);
    epsilon = L(1,1);
    % Normalize wavefunction
    norm = 1/sqrt(h);
    f = K(:,1) * norm * sign(K(1,1));
    f = f .* 1./(rValues * sqrt(4*pi));
    
    % Compute energy
    ns = f.^2;
    E = 2*epsilon - sum(4*pi*Vh.*rValues.^2 .* ns * h);
    
    energy = [energy, E];
    disp(['Iteration: ' num2str(i) ', E = ' num2str(E)]);
    i=i+1;
end
% Plot results
figure(1)
clf;
hold on;
plot(rValues,4*pi*rValues.^2 .*ns)
plot(rValues, 4*rValues.^2 * 8 .* exp(-4*rValues),'--')
plot(rValues, 4*rValues.^2 * (27/16)^3 .* exp(-2*(27/16)*rValues),'.--')
hold off;
legend('Simulation','Theoretical, Z=2','Theoretical Z=27/16')
xlabel('r - distance from nucles [a_0]')
ylabel('p(r) - probability of findning electron at r')
title(['Helium, simulated ground state energy = ' num2str(E) ' H.E'])