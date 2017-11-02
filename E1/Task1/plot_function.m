% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data1 = importdata('function_phi0.dat');
data2 = importdata('function_phipihalf.dat');

%plot
figure;
hold on;
plot(data1(:,1),data1(:,2),'-');
plot(data2(:,1),data2(:,2),'-.');
legend(["\phi = 0", "\phi = \pi/2"]);

% labels
xlabel('$t$','Interpreter','LaTeX');
ylabel('$h(t)$','Interpreter','LaTeX');



