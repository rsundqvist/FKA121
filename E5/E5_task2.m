clc, clear all, close all;
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

% Set up Schrödinger equation to compute energy
potential = zeros(nbrOfPoints);
for i=2:nbrOfPoints-1
   potential(i,i) = 1/rValues(i); 
end
% Hamiltonian
H = -0.5*A - potential;
% Compute energy and new wave function
[K,L] = eig(H);
groundStateEnergy = L(1,1)

% Normalize wavefunction
norm = 1/sqrt(h * sum(K(:,1).^2));
K(:,1) = K(:,1) * norm;

%% Plot the simulated wavefunction
figure(1)
clf;
hold on;
plot(rValues,f(rValues))
plot(rValues,-K(:,1),'--','LineWidth',4)
hold off;
legend('Theoretical','Simulated')
title(['Ground state Hydrogen, simulated energy = ' num2str(groundStateEnergy)])
xlabel('Distance from nuclues [a_0]')
ylabel('Radial wavefunction')