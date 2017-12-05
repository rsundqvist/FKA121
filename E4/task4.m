close all;
clear all;
clc;

N = 4;
Q = [];
V = [];
for i = 0:N-1
    v = load(['vel_tau' num2str(i) '.dat']);
    V = [V v];
    q = load(['offset_tau' num2str(i) '.dat']);
    Q = [Q q];
end

figure(1);

subplot(2, 1, 2);
hold on;
for i = 1:N
    histogram(V(:,i), 'normalization', 'probability');
end
title('Velocity distribution');
xlabel('V/V_{th}');
ylabel('Probability');
legend({'$0.5\tau$', '$\tau$', '$2\tau$', '$10\tau$'}, 'interpreter', 'latex');


subplot(2, 1, 1);
hold on;
for i = 1:N
    histogram(Q(:,i), 'normalization', 'probability');
end
title('Offset distribution');
xlabel('Q [m]');
ylabel('Probability');
legend({'$0.5\tau$', '$\tau$', '$2\tau$', '$10\tau$'}, 'interpreter', 'latex');
