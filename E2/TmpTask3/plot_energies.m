<<<<<<< HEAD
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
plot(data(:,1),data(:,4))
plot(data(:,1),data(:,5))
plot(data(:,1),data(:,6))
plot(data(:,1),data(:,7))
hold off
=======
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
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');
<<<<<<< HEAD
%axis([0,250,0,35])
=======
axis([0 250 0 35])
>>>>>>> 6613cc9cc7b5a7b9919a5e3174fc5ceed1602f32
