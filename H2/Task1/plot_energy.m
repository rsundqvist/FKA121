data10 = load('energy_stats10.dat');
data90 = load('energy_stats90.dat');
data2 = load('energy_traj.dat');


figure

subplot(2, 1, 1);
mu = data10(:,1);
sigma = sqrt(data10(:,2));
hold on;
grid on
plot(mu, 'LineWidth', 2);
plot(mu+sigma, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
plot(mu-sigma, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
title('Energy, accept rate < 10%');

subplot(2, 1, 2);
mu = data90(:,1);
sigma = sqrt(data90(:,2));
hold on;
grid on
plot(mu, 'LineWidth', 2);
plot(mu+sigma, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
plot(mu-sigma, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
title('Energy, accept rate > 90%');


for i = 1:3:10
    %plot(data2(i, :), 'LineWidth', 1, 'LineStyle', '-.');
end

legend('\mu_{E_L}', '\mu_{E_L} \pm \sigma');