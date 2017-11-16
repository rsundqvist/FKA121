close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;
plot(t,data(:,2))
plot(t,data(:,3))
plot(t,data(:,4))
hold off;
legend('Potential energy','Kinetic energy','Total energy')
