%% Plot vol
close all;
clear all;

data = load('velcorr_973.dat');
%data2 = load('msd_773.dat');



t = data(1:end,1);
figure(1);


hold on;
plot(t(1:100),data(1:100,14))
%plot(t,data2(1:end,14))

xlabel('Time [ps]')
ylabel('msd')
%legend('\Tau_{eq} = 973K', '\Tau_{eq} = 773K');

sz = 60000;
W = linspace(0, 1000,sz);
phi = zeros(1, sz);
for i=1:sz
    w = W(i);
    tmpData = data(1:sz,14).* cos(w*t(1:sz));
    phi(i) = 2 * trapz(tmpData);
end
plot(t(1:sz),phi)
%Q = 0.33333333*trapz(data(1:100,14))

hold off;