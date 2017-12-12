% Initialize grid
nbrOfPoints = 10000;
h = 0.01; % grid spacing
r0 = 0;
rValues = r0 + (1:nbrOfPoints)*h; 

% Set electron density according to the ground state of hydrogen

ns = 4 * exp(-2*rValues);

% Initialize matrix
A = zeros(nbrOfPoints);
for i=2:nbrOfPoints-1 % Inner points
    A(i,i+1) = 1/(h^2);
    A(i,i-1) = 1/(h^2);
    A(i,i) = - 2/(h^2);
end

% Boundary conditions
A(1,1) = -2/(h^2);
A(1,2) = 1/(h^2);
A(end,end) = -2/(h^2);
A(end,end-1) = 1/(h^2);

c = zeros(nbrOfPoints,1);
c(end) = 1;
% Set matrix equation Ax = b + c
U = A\ (ns' + c);

% Solve Hartree potential and compare to analytical solution
Vh = @(r) 1./r - (1 + 1./r).*exp(-2*r);
figure(1)
hold on
plot(Vh(rValues),rValues)
plot(rValues,U)