zclear all, close all;
% Initialize grid
nbrOfPoints = 10000;
r0 = 0.001;
rEnd = 10;
rValues = linspace(r0,rEnd,nbrOfPoints); 
h = rValues(2) - rValues(1); % grid spacing

% Set electron density according to the ground state of hydrogen
electronDensity = @(r) exp(-2.*r)/pi;
% Multiply ns with r, since it will be needed in the equation later
ns = rValues .* electronDensity(rValues);%4 * exp(-2*rValues);

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

% Set matrix equation Ax = b + c
U = A\(- 4*pi*ns' + c);
Vh_sim = U./rValues';

% Solve Hartree potential and compare to analytical solution
Vh = @(r) 1./r - (1 + 1./r).*exp(-2*r);
figure(1)
hold on
plot(rValues, Vh(rValues))
plot(rValues, Vh_sim,'--','Linewidth',3)

legend('Hartree Potential','Simulated Potential')
title('Electronic potential for Hydrogen')
xlabel('Distance to nuclues [a_0]')
ylabel('Hartree energy')
