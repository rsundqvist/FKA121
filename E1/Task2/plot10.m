
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = 1, \phi = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data file
dataps = importdata('dat/powerspectrum_1_0.dat');
dataf = importdata('dat/function_1_0.dat');
%plot
pos = pos+1;
subplot(rows,cols,pos);
plot(dataps(:,1),dataps(:,2),'-');
title('f = 1, \phi = 0');
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