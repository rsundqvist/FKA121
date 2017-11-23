%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;

P_mean = mean(data(600:end,5))
plot(t,data(:,5))
plot([0, max(t)],[P_mean P_mean], 'k');

T_mean = mean(data(600:end,6))
plot(t,data(:,6))
plot([0, max(t)],[T_mean T_mean], 'k');


legend('T [K]', '<P>_t', 'P [eV/Å^3]', '<P>_t');
xlabel('Time [ps]')
hold off;