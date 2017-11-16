%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;
plot(t,data(:,1), '-.')
plot(t,data(:,2), '.')
plot(t,data(:,3))
hold off;

legend('E_p', 'E_k', 'E_{tot}');