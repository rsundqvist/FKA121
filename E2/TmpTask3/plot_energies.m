% load the data file
data = importdata('disp.dat');

%plot
figure;
plot(data(:,1),data(:,2),data(:,1),data(:,3),data(:,1),data(:,4),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Energy / [dim. unit]');



