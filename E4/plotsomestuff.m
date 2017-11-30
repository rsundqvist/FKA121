close all;
clear all;
clc;
data = load('e4.dat');

lim = 1:size(data, 1);
t = data(lim, 1);
mu_q = data(lim, 2);
sigma_q = sqrt(data(lim, 3));
mu_v = data(lim, 4);
sigma_v = sqrt(data(lim, 5));

traj_q = data(lim, 6:10);
vel_q = data(lim, 11:15);

figure;
subplot(2,1,1);
hold on;
plot(t, mu_q);
plot(t, mu_q + sigma_q, 'k');
plot(t, mu_q - sigma_q, 'k');
legend('\mu','\mu \pm \sigma');
title('Average offset');
xlabel('Time [\mus]');

subplot(2,1,2);
hold on;
plot(t, mu_v);
plot(t, mu_v + sigma_v, 'k');
plot(t, mu_v - sigma_v, 'k');
legend('\mu','\mu \pm \sigma');
title('Average speed');
xlabel('Time [\mus]');

figure;
subplot(2,1,1);
hold on;
plot(t, traj_q(:,1));
plot(t, traj_q(:,2));
plot(t, traj_q(:,3));
plot(t, traj_q(:,4));
plot(t, traj_q(:,5));
ylabel('q [\mum]');
xlabel('Time [\mus]');
subplot(2,1,2);
hold on;
plot(t, vel_q(:,1));
plot(t, vel_q(:,2));
plot(t, vel_q(:,3));
plot(t, vel_q(:,4));
plot(t, vel_q(:,5));
ylabel('v [\mum]');
xlabel('Time [\mus]');