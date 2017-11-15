% plot the energies and displacements
close all;
clear all;
% load the data file
data1 = importdata('displacements.dat');
data2 = importdata('energy.dat');

%plot displacements
data = data1;
figure(1);
hold on
plot(data(:,1),data(:,2))
plot(data(:,1),data(:,3))
plot(data(:,1),data(:,4))
plot(data(:,1),data(:,5))
plot(data(:,1),data(:,6))
hold off

% labels
xlabel('Time / [dim. unit]');
ylabel('Displacements / [dim. unit]');


%plot energies
data = data2;
figure(2);
hold on
plot(data(:,1),data(:,2))
plot(data(:,1),data(:,3))
hold off

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');
axis([0,250,0,35])