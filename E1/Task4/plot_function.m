% plot the powerspectrum
% Created by Martin Gren 2014-10-25.
close all;

% load the data file
data1 = importdata('function.dat');
data2 = importdata('powerspectrum.dat');

%plot
figure;
hold on;
<<<<<<< HEAD
subplot(2,1,1)
plot(data1(:,1),data1(:,2),'-');
=======
subplot(2, 1, 1);
plot(data1(:,1),data1(:,2),'-');
legend('function')
subplot(2, 1, 2);
plot(data2(:,1),data2(:,2),'-.');
legend('powerspec');

% labels
>>>>>>> 42c4993f498120e886ce1782b145921a49803590
xlabel('$t$','Interpreter','LaTeX');
ylabel('$h(t)$','Interpreter','LaTeX');
subplot(2,1,2)
plot(data2(:,1),data2(:,2),'-.');
xlabel('$f$','Interpreter','LaTeX');
ylabel('$Pn$','Interpreter','LaTeX');
%legend('function', 'powahspec');
