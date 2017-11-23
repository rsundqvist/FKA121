%% Plot energies
close all;
clear all;

data = load('simulation2.dat');

t = data(:,1);
equibEnd = find(t>600, 1);
figure(1);
hold on;

X1 = data(:,7:9);
X2 = data(:,10:12);

plot3(X1(equibEnd:end,1), X1(equibEnd:end,2), X1(equibEnd:end,3));
plot3(X2(equibEnd:end,1), X2(equibEnd:end,2), X2(equibEnd:end,3));

legend({'$\vec{r}_1$', '$\vec{r}_{70}$'},'Interpreter','latex');
grid on
hold off;