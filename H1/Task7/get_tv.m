data = load('velcorr_773.dat');

sz = 50000;
t = data(1:sz,1);
mean_vel = data(1:sz,14);

plot(t(1:100), mean_vel(1:100));
dlmwrite('tv_773.dat', [t mean_vel]);