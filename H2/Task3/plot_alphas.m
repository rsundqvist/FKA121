data = load('alphas.dat');
%%
clf;
hold on
plot(data(:,1),data(:,2))
plot(data(:,1),data(:,2)+data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(data(:,1),data(:,2)-data(:,3), 'Color', 'k', 'LineStyle', '--')
hold off

xlabel('\alpha')
ylabel('Ground state energy E_0 [Hartree energy]')
legend('<E_L>','<E_L> \pm \sigma')

I = find(min(data(:,2))==data(:,2), 1)
data(I,1)