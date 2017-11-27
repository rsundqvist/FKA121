%% Plot vol
close all;
clear all;

data = load('msd_973.dat');
data2 = load('msd_773.dat');

t = data(1:end,1);
figure(1);


hold on;
plot(t,data(1:end,14))
plot(t,data2(1:end,14))

xlabel('Time [ps]')
ylabel('msd')
%legend('\Tau_{eq} = 973K', '\Tau_{eq} = 773K');
legend('973', '773');

hold off;