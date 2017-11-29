
clear all;
<<<<<<< HEAD
%data = load('allpos_solid.dat');
%pos = data(:, 2:end);
data = load('allpos_liquid.dat');
pos = data2(:, 2:end);
=======
disp('Loading solid data.');
data = load('allpos_solid.dat');
pos = data(:, 2:end);
disp('Loading liquid data.');
data2 = load('allpos_liquid.dat');
pos2 = data2(:, 2:end);
disp('Data loaded.');
>>>>>>> 5d8b81613866e3c6e7c397a86db10e2d90585cbc
