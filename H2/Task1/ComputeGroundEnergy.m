clc, close all, clear all;
data = load('trajectory.dat');

alpha = 0.1;

X1 = data(:,1);
Y1 = data(:,2);
Z1 = data(:,3);

X2 = data(:,4);
Y2 = data(:,5);
Z2 = data(:,6);

energySum = 0;
energy = zeros(1,length(data));

for i=1:length(data)
    r1 = [X1(i),Y1(i),Z1(i)];
    r2 = [X2(i),Y2(i),Z2(i)];
    energy(i) = computeLocalEnergy(r1,r2,alpha);
    energySum = energySum + computeLocalEnergy(r1,r2,alpha);
end

energySum/length(data)
figure
plot(energy)
title('energy')

%% Sample from markov chain
S = 37;
sampleEnergy = energy(1:S:end); 
figure
plot(sampleEnergy)
title('sampleEnergy')

%% Plot std over iterations
figure
data = sampleEnergy;
stds = zeros(1,length(data));
for i=1:length(data)
   stds(i) = std(energy(1:i)); 
end

plot(stds)
title('stds')