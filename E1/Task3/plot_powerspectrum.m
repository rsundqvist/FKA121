% plot the powerspectrum
% Created by Martin Gren 2014-10-25.
close; % close figs
% load the data file
data = importdata('powerspectrum258.dat');

%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Frequency','Interpreter','LaTeX');
ylabel('Powers pectrum','Interpreter','LaTeX');
title('f = 2, \psi = 0, n = 258');


