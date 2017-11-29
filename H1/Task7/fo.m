data = load('velcorr_973.dat');

sz = 50000;
t = data(1:sz,1);
mean_vel = data(1:sz,14);

dlmwrite('tv_973.dat', [t mean_vel]);