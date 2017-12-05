close all;
clear all;

data = load('block_average.dat');

t = 15000;
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));
plot([0 lim], [10 10]);

legend('\Phi_k', 's = 15');