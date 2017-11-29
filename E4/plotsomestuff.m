data = load('e4.dat');

lim = 10;
t = data(1:lim, 1);
mu_q = data(1:lim, 2);
sigma_q = sqrt(data(1:lim, 3));
mu_v = data(1:lim, 4);
sigma_v = sqrt(data(1:lim, 5));

traj_q = data(1:lim, 6:end);

figure;
subplot(2,1,1);
hold on;
plot(t, mu_q);
plot(t, mu_q + sigma_q, ':');
plot(t, mu_q - sigma_q, ':');
legend('foo','bar','baz');