%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;

X1 = data(:,7:9);
X2 = data(:,10:12);

plot3(X1(100:end,1), X1(100:end,2), X1(100:end,3));
plot3(X2(100:end,1), X2(100:end,2), X2(100:end,3));
hold off;