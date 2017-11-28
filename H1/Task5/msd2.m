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
msd = computeMSD(pos,maxStep,N);

plot(msd);

disp('Beep!');
beep;