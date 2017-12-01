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
figure(2)
hold on
plot(r1)
plot(r2)
plot(r);
hold off

figure(3)
%r = [r1; r2];
histogram(r,'Normalization','pdf');

figure(4);
theta = data(:,7); %angle between electrons
histogram(theta,'Normalization','pdf');