% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('powerspectrum.dat');

%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Frequency / [dim. unit]','Interpreter','LaTeX');
ylabel('Powerspectrum / [dim. unit]','Interpreter','LaTeX');

% axis limits
xlim([-200,200]);

%% Frequencies -> cm^-1
c = 2.99792458*10^8; % speed of light

f1 = 38.91*10^12;
f2 = 74.77*10^12;

k1 = f1/c * 10^(-2)
k2 = f2/c * 10^(-2)

%% Task 8 calculate stuff

k = 1600; %N/m
m = 16 * 1.66053904*10^(-27);
M = 12 * 1.66053904*10^(-27);

lamda1 = 0
lambda2 = -k*(-2*m - M)/(m*M)
lambda3 = k/m

omega1 = 2*sqrt(k/M)*sin(k*pi/(2*(3)));
omega2 = 2*sqrt(k/m)*sin(k*pi/(2*(3)));
f1a = omega1/(2*pi);
f2a = omega2/(2*pi);

k1a = f1a/c * 10^(-2)
k2a = f2a/c * 10^(-2)
