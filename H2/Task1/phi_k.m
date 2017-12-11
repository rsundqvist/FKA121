close all;
clear all;

data = load('phi_k.dat');

t = 30;
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));
plot([0 lim], [exp(-2) exp(-2)]);

legend('\Phi_k', 'e^{-2}');
grid on;
grid minor

tmp = find(data(:,2)<exp(-2), 1) - 1;
plot(tmp,data(tmp+1,2),'o')