%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;
plot(t,data(:,2), '-.')
plot(t,data(:,3), '.')
plot(t,data(:,4))
hold off;

legend('E_p', 'E_k', 'E_{tot}');

%% Compute temperature
N = 256;
kB = 0.0000080000617330;
Ek = data(1:end,3);

T = 2/(kB*3*N) * mean(Ek)

%%
figure(2)
plot(data(:,1),data(:,5))

figure(3)
plot(data(:,1),data(:,6))