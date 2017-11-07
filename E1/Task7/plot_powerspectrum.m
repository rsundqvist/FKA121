% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('powerspectrum.dat');

%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Frequency / [dim. unit]','Interpreter','LaTeX');
ylabel('Powerspectrum / [dim. unit]','Interpreter','LaTeX');

% axis limits
xlim([-200,200]);

%%
c = 2.99792458*10^8; % speed of light

f1 = 38.91*10^12;
f2 = 74.77*10^12;

k1 = f1/c * 10^(-2)
k2 = f2/c * 10^(-2)