close all;
clear all;

data = load('block_average.dat');

s = 10;
t = 2000;
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));
plot([0 t], [s s]);

legend('\Phi_k', ['s = ' num2str(s)]);