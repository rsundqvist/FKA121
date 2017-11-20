close all;
clear all;
data = load('hist.dat');

figure(1)
hold on
histogram(data,'Normalization','pdf')
n = length(data);
title(['Number of points N = ' num2str(n)])
% Plot distribution
figure(2)
p = @(x) pi/2 * sin(pi*x); 
X = linspace(0,1,1000);
plot(X,p(X))
