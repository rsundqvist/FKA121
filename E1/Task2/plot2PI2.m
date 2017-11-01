


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = 2, \phi = \pi/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data file
dataps = importdata('dat/powerspectrum_2_PI2.dat');
dataf = importdata('dat/function_2_PI2.dat');
%plot
pos = pos+1;
subplot(rows,cols,pos);
plot(dataps(:,1),dataps(:,2),'-');
title('f = 2, \phi = \pi/2');
% labels
xlabel('Frequency','Interpreter','LaTeX');
ylabel('Powerspectrum','Interpreter','LaTeX');

%plot
pos = pos+1;
subplot(rows,cols,pos);
plot(dataf(:,1),dataf(:,2),'-');
% labels
xlabel('t','Interpreter','LaTeX');
ylabel('h(t)','Interpreter','LaTeX');