% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data1 = importdata('function_phi0_f2.dat');
data2 = importdata('function_phiPiHalf_f1.dat');
data3 = importdata('function_phi0_f1.dat');
%plot
figure;
hold on;
plot(data1(:,1),data1(:,2),'-');
plot(data2(:,1),data2(:,2),'-.');
plot(data3(:,1),data3(:,2),'-.');
legend('\phi = 0, f = 2', '\phi = \pi/2, f = 1', '\phi = 0, f = 1');

% labels
xlabel('$t$','Interpreter','LaTeX');
ylabel('$h(t)$','Interpreter','LaTeX');



