%% Plot energies
close all;
clear all;

data = load('msd_973.dat');

t = data(:,1);
figure(1);

subplot(2,1,1);
hold on;
equibEnd = find(t>600, 1)
T_mean = mean(data(equibEnd:end,5))
plot(t,data(:,5), 'r')
plot([0, max(t)],[T_mean T_mean], 'k');
legend('T [K]', '<T>');
xlabel('Time [ps]')
ylabel('Temperature [K]')

subplot(2,1,2);
hold on;
P_mean = mean(data(equibEnd:end,6))
plot(t,data(:,6))
plot([0, max(t)],[P_mean P_mean], 'k');
plot([0, max(t)],[0.000000633 0.000000633], ':');
legend('P [eV/Å^3]', '<P>', 'P_{eq}');
xlabel('Time [ps]')
ylabel('Pressure [eV/�^3]')


hold off;