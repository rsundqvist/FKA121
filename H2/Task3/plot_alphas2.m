data = load('alphas.dat');

N = 10^5;
alphas = data(:,1);
energy = data(:,2);
err = data(:,2)/sqrt(N);

figure;
plot(alphas, energy);
ylabel('E [E_{hartree}]');
xlabel('\alpha');
%errorbar(alphas, energy);

Emin = min(energy(:));
I = find(Emin==energy(:), 1);
amin = alphas(I)
title(['E_{min} = ' num2str(Emin) ', \alpha_{min} = ' num2str(amin)]);