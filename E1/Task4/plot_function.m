% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data1 = importdata('function.dat');
data2 = importdata('powerspectrum.dat');

%plot
figure;
hold on;
plot(data1(:,1),data1(:,2),'-');
plot(data2(:,1),data2(:,2),'-.');
legend('function', 'powahspec');

% labels
xlabel('$t$','Interpreter','LaTeX');
ylabel('$h(t)$','Interpreter','LaTeX');



