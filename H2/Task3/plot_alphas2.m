data = load('alphas.dat');

N = 10^5;
alphas = data(:,1);
energy = data(:,2);
err = ones(size(alphas)) * std(alphas)/sqrt(N);

avg = mean(energy)
dev = std(energy)

figure;
hold on
grid on
p = polyfit(alphas,energy,2);
yp = polyval(p, alphas, '-.');
errorbar(alphas, energy, err);
%plot(alphas, energy);
plot(alphas, yp);
ylabel('E(\alpha) [E_{hartree}]');
xlabel('\alpha');

Emin = min(energy(:))
I = find(Emin==energy(:), 1);
amin = alphas(I)

I = find(min(yp)==yp(:), 1);
min(yp)
amin = alphas(I)
%title(['E_{min} = ' num2str(Emin) ', \alpha_{min} = ' num2str(amin)]);
title('Energy values');
legend('Simulation energy', '2nd-degree polynomial fit');