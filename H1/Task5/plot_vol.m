%% Plot vol
close all;
clear all;

data = load('vol.dat');

t = data(1:end,1);
figure(1);


hold on;
plot(t,data(1:end,13))

xlabel('Time [ps]')
ylabel('Supercell length [Å]')

hold off;