% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('powerspectrumPI2.dat');

%plot
figure;
subplot(2,1,1);
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Frequency','Interpreter','LaTeX');
ylabel('Powerspectrum','Interpreter','LaTeX');

% load the data file
data = importdata('functionPI2.dat');

%plot
subplot(2,1,2);
plot(data(:,1),data(:,2),'-');

% labels
xlabel('t','Interpreter','LaTeX');
ylabel('h(t)','Interpreter','LaTeX');



