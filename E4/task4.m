close all;
clear all;
clc;

Theo = 297;
kB = 1.38065e-23;
m = 30.131e-15;
v0 = 2.0e-3;
v_th = 3.689039e-04;

Tau_a = 48.5e-6;  % relaxation time in microseconds, case A
Tau_b = 147.3e-6; % case B
Tau = Tau_a;
eta = 1/Tau;
    
t = [0.1*Tau, 0.5*Tau, Tau, 2*Tau, 10*Tau];
vRange = linspace(-8, 8, 1000)*v_th;
f = @(v,t) sqrt( m/(2*pi*kB*Theo*(1-exp(-2*eta*t))) ) .* exp( -m*(v-v0*exp(-eta*t)).^2 ./ (2*kB*Theo*(1-exp(-2*eta*t))) );

N = 5;
Q = [];
V = [];
Theo = []; %Theoretical
for i = 0:N-1
    v = load(['vel_tau' num2str(i) '.dat']);
    V = [V v];
    q = load(['offset_tau' num2str(i) '.dat']);
    Q = [Q q];
end

figure(1);

subplot(3, 1, 2);
hold on;
for i = 1:N
    histogram(V(:,i), 'normalization', 'probability');
end

subplot(3, 1, 3);
hold on
for i = 1:N
    plot(vRange/v_th, f(vRange, t(i)), 'Color', 'k', 'LineStyle', ':');
end
title('Velocity distribution');
xlabel('V/V_{th}');
ylabel('Probability');
legend({'$0.1\tau$', '$0.5\tau$', '$\tau$', '$2\tau$', '$10\tau$', 'Theoretical'}, 'interpreter', 'latex');


subplot(3, 1, 1);
hold on;
for i = 1:N
    histogram(Q(:,i), 'normalization', 'probability');
end

title('Offset distribution');
xlabel('Q [m]');
ylabel('Probability');
legend({'$0.1\tau$', '$0.5\tau$', '$\tau$', '$2\tau$', '$10\tau$'}, 'interpreter', 'latex');
