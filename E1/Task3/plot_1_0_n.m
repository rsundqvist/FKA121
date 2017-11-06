%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = 1, \phi = 0, n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data file
dataps = importdata(sprintf('dat/powerspectrum_f1_p0_n%d.dat', n));
dataf = importdata(sprintf('dat/function_f1_p0_n%d.dat', n));
pos = pos+1;
subplot(rows,cols,pos);
plot(dataps(:,1),dataps(:,2),'-');
title(sprintf('f = 1, \\phi = 0, n = %d',n));
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