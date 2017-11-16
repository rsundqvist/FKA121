close all;
clear all;

% load the data file
data1 = importdata('disp.dat');
data2 = importdata('energy.dat');
%plot
data = data1;
figure(1);
plot(data(:,1),data(:,2),data(:,1),data(:,3),data(:,1),data(:,4),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');

% Plot energy
data = data2;
figure(2);
plot(data(:,1),data(:,2),data(:,1),data(:,3),data(:,1),data(:,4),data(:,1),data(:,5),data(:,1),data(:,6),data(:,1),data(:,7),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');
axis([0 250 0 35])