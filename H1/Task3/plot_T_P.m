%% Plot energies
close all;
clear all;

data = load('energy.dat');

t = data(:,1);
figure(1);
hold on;

equibEnd = find(t>600, 1)
P_mean = mean(data(equibEnd:end,5))
plot(t,data(:,5))
plot([0, max(t)],[P_mean P_mean], 'k');

T_mean = mean(data(equibEnd:end,6))
plot(t,data(:,6))
plot([0, max(t)],[T_mean T_mean], 'k');


legend('T [K]', '<P>', 'P [eV/Å^3]', '<P>');
xlabel('Time [ps]')
hold off;