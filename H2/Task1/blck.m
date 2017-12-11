close all;
clear all;

data = load('block_average.dat');

s = 10;
t = size(data,1);
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));
plot([0 lim], [12 12]);

xlabel('Block size')
ylabel('Block average')
legend('\Phi_k', 's = 12');
plot([0 t], [s s]);

legend('\Phi_k', ['s = ' num2str(s)]);
