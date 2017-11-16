%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;
plot(t,data(:,4))
hold off;

%% Plot trajectory

figure(2)
plot3(data(:,),Y1,Z1)