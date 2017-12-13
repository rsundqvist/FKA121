clear all, close all;
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
for i=3:nbrOfPoints-2 % Inner points
    A(i,i+1) = 1/(h^2);
    A(i,i-1) = 1/(h^2);
    A(i,i) = - 2/(h^2);
end

% Boundary conditions
%A(1,1) = -2/(h^2);
%A(1,2) = 1/(h^2);
%A(end,end) = -2/(h^2);
%A(end,end-1) = 1/(h^2);
A(1,2) = 1/(h^2);
A(2,2) = - 2/(h^2);
A(2,3) = 1/(h^2);
A(end,end-1) = 1/(h^2);
A(end-1, end-1) = - 2/(h^2);
A(end-1, end-2) = 1/(h^2);

c = zeros(nbrOfPoints,1);
c(end) = - 2/(h^2);
c(end-1) = 1/(h^2);

% Set matrix equation Ax = b + c
U = pinv(A)*(- 4*pi*ns' + c);

% Solve Hartree potential and compare to analytical solution
Vh = @(r) 1./r - (1 + 1./r).*exp(-2*r);
figure(1)
hold on
plot(rValues, Vh(rValues))
plot(rValues,U./rValues')