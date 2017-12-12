close all;
clear all;

data = load('block_average.dat');

s = round(mean(data(250:end, 2)))
t = size(data,1);
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));

xlabel('Block size')
ylabel('Block average')
plot([2 t+2], [s s], 'Color', 'k', 'LineStyle','-.');

legend('\Phi_k', ['s = ' num2str(s)]);

title('Box average');
