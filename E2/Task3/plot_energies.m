%% Plot the energies in 5 particles
clear all;
close all;

% load the data file
data = importdata('energy.dat');

figure;
hold on
plot(data(:,1),data(:,2));
plot(data(:,1),data(:,3));
plot(data(:,1),data(:,4));
plot(data(:,1),data(:,5));
plot(data(:,1),data(:,6));
plot(data(:,1),data(:,7), '-.');
hold off
legend('E_1(t)','E_2(t)','E_3(t)','E_4(t)','E_5(t)', 'E_{tot}(t)')
xlabel('Time')
ylabel('Energy')
axis([0,250,0,40])