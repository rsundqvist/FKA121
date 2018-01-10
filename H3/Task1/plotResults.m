data = load('simulation.dat');
%posData = load('position.dat');

nWalkers = data(:,2:2:end);
energy = data(:,1:2:end-1);

meanEnergy = mean(energy);
meanWalkers = mean(nWalkers);
stdEnergy = std(energy);
stdWalkers = std(nWalkers);

figure(1)
clf;
hold on;
plot(meanEnergy)
plot(meanEnergy+stdEnergy,'red')
plot(meanEnergy-stdEnergy,'red')
hold off
figure(2)
clf;
hold on;
plot(meanWalkers)
plot(meanWalkers+stdWalkers,'red')
plot(meanWalkers-stdWalkers,'red')
hold off;

figure(3)
clf;
hold on
plot(nWalkers(ceil(rand*size(nWalkers,1)),:))
plot(nWalkers(ceil(rand*size(nWalkers,1)),:))
plot(nWalkers(ceil(rand*size(nWalkers,1)),:))
plot(nWalkers(ceil(rand*size(nWalkers,1)),:))
plot(nWalkers(ceil(rand*size(nWalkers,1)),:))
plot(nWalkers(2,:))
plot(nWalkers(4,:))
plot(nWalkers(88,:))
hold off

figure(4)
clf;
hold on;
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
plot(energy(ceil(rand*size(energy,1)),:))
hold off;

%tmpArray = [];
%for i=1:size(posData,1)*size(posData,2)
%   if isnan(posData(i)) == 0
%      tmpArray = [tmpArray; posData(i)]; 
%   end
%end
%posData = tmpArray;

%figure(5)
%plot(mean(posData))
