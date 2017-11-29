%% Plot vol
close all;

%loadAllpos;

dt = 1;
timeSteps = size(data, 1);

beep on;

%h = waitbar(0, 'Progress - overall');
%h2 = waitbar(0, 'Progress - particle');
maxStepSize = 1; % Step size if amplified by ir in MD_main.c!
N = 256;
maxStep = size(pos, 1);

disp('begin computeMSD');
%msd_solid = computeMSD(pos,maxStep,N);
%msd_liquid = computeMSD(pos2,maxStep,N);

figure;
hold on;
%plot(msd_solid, 'DisplayName', 'Solid');
%plot(msd_liquid, 'DisplayName', 'Liquid');

t = 1:2000;
polyfit(t, msd_liquid(t), 1)
D = 1./(6*t) .* msd_liquid(t);

plot(t, D);
disp('Beep!');
beep;


