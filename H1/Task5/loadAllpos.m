
clear all;
disp('Loading solid data.');
data = load('allpos_solid.dat');
pos = data(:, 2:end);
disp('Loading liquid data.');
data2 = load('allpos_liquid.dat');
pos2 = data2(:, 2:end);
disp('Data loaded.');