clc, close all, clear all;
data = load('trajectory.dat');

figure(1);
hold on;
grid on;
plot3(data(:,1),data(:,2),data(:,3))
plot3(data(:,4),data(:,5),data(:,6))
[x, y, z] = sphere(256);
h = surfl(x, y, z); 
set(h, 'FaceAlpha', 0.7)
shading interp
axis([-3 3 -3 3 -3 3]);

r1 = sqrt(data(:,1).^2 + data(:,2).^2 + data(:,3).^2);
r2 = sqrt(data(:,4).^2 + data(:,5).^2 + data(:,6).^2);
figure(2)
hold on
plot(r1)
plot(r2)
hold off

figure(3)
r = [r1; r2];
histogram(r,'Normalization','pdf');