clc, close all, clear all;
data = load('trajectory.dat');

figure(1);
hold on;
grid on;

X1 = data(:,1);
Y1 = data(:,2);
Z1 = data(:,3);

X2 = data(:,4);
Y2 = data(:,5);
Z2 = data(:,6);

Xm = (X1+X2)/2;
Ym = (Y1+Y2)/2;
Zm = (Z1+Z2)/2;

plot3(X1, Y1, Z1);
plot3(X2, Y2, Z2);
plot3(Xm, Ym, Zm);

[x, y, z] = sphere(256);
h = surfl(x, y, z); 
set(h, 'FaceAlpha', 0.7)
shading interp
axis([-3 3 -3 3 -3 3]);

r1 = sqrt(data(:,1).^2 + data(:,2).^2 + data(:,3).^2);
r2 = sqrt(data(:,4).^2 + data(:,5).^2 + data(:,6).^2);
r = (r1+r2)/2;

figure
%subplot(2,1,1)
yyaxis left
h = histogram(r,'Normalization','pdf');
xlabel('r [a_0]')
ylabel('Probability distribution')
title('Variational Monte Carlo')
axis([0,3,0,1.2])

%subplot(2,1,2)
yyaxis right
hold on
rtmp = linspace(0,5.5,1000);
f1 = @(x) 2^3 * 4 * x.^2 .* exp(-2*2*x);
f2 = @(x) (27/16)^3 * 4 * x.^2 .* exp(-2*(27/16)*x);
plot(rtmp,f1(rtmp), 'Color', 'k', 'LineStyle', '-')
plot(rtmp,f2(rtmp), 'Color', 'k', 'LineStyle', '--')
legend({'Simulation', '$p(r), Z = 2$', '$p(r), Z = \frac{27}{16}$'}, 'Interpreter', 'latex')
xlabel('r [a_0]')
ylabel('Probability')
title('Central-field approximation')
axis([0,3,0,1.2])
hold off
%% Theta
figure(4);
theta = data(:,7); %angle between electrons
x = cos(theta);
histogram(theta,'Normalization','pdf');
figure(5)
plot(x)