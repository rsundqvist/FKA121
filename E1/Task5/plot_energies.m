% plot the energies
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('energy.dat');

%plot
figure;
plot(data(:,1),data(:,2),data(:,1),data(:,3),data(:,1),data(:,4),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');

% legend
legend('Total energy','Potential energy','Kinetic energy');

% axis limits
xlim([0,50]);
