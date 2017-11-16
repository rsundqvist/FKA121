%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;
plot(t,data(:,5))
plot(t,data(:,6))
hold off;