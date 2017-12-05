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
energia = zeros(1,length(data));

for i=1:length(data)
    r1 = [X1(i),Y1(i),Z1(i)];
    r2 = [X2(i),Y2(i),Z2(i)];
    energia(i) = computeLocalEnergy(r1,r2,alpha);
    energySum = energySum + computeLocalEnergy(r1,r2,alpha);
end

energySum/length(data)
plot(energia)
%% Sample from markov chain
S = 37;
sampleEnergia = zeros(1,floor(length(energia)/S));
for i=1:floor(length(energia)/S)
   sampleEnergia(i) = energia(i*S); 
end

plot(sampleEnergia)

%% Plot std over iterations
data = sampleEnergia
stds = zeros(1,length(data));
for i=1:length(data)
   stds(i) = std(energia(1:i)); 
end

plot(stds)