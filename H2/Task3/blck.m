close all;
clear all;

data = load('block_average.dat');

s = 10;
t = 2000;
hold on;
plot(data(1:t,1), data(1:t,2));
lim = max(data(1:t,1));
<<<<<<< HEAD
plot([0 lim], [12 12]);

xlabel('Block size')
ylabel('Block average')
legend('\Phi_k', 's = 12');
=======
plot([0 t], [s s]);

legend('\Phi_k', ['s = ' num2str(s)]);
>>>>>>> aa603691b4d55d3542dda8cf051ddcff1f4b469e
