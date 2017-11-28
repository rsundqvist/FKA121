%% Plot vol
close all;

%loadAllpos;

dt = 1;
timeSteps = size(data, 1);

beep on;

h = waitbar(0, 'Progress - overall');
h2 = waitbar(0, 'Progress - particle');
maxStepSize = 1; % Step size if amplified by ir in MD_main.c!

ir = 5;
N = 256;

msd_solid = zeros(maxStepSize, timeSteps);
for stepSize=1:maxStepSize
    waitbar(0.5*stepSize/maxStepSize, h);
    for i=1:N
        waitbar(i/N, h2);
        for j = 1:timeSteps-stepSize
            i1 = 1+(i-1)*3; i2 = 3+(i-1)*3; % x to z for the particle
            r1 = pos(j+stepSize,i1:i2); % Position after stepSize iterations
            r0 = pos(j,i1:i2); % Base position
            d = r1-r0;
            msd_solid(stepSize, j) = msd_solid(stepSize, j) + dot(d, d); % Add displacement of the particle
        end
    end
end

msd_liquid = zeros(maxStepSize, timeSteps);
for stepSize=1:maxStepSize
    waitbar(0.5 + stepSize/maxStepSize, h);
    for i=1:N
        waitbar(i/N, h2);
        for j = 1:timeSteps-stepSize
            i1 = 1+(i-1)*3; i2 = 3+(i-1)*3; % x to z for the particle
            r1 = pos2(j+stepSize,i1:i2); % Position after stepSize iterations
            r0 = pos2(j,i1:i2); % Base position
            d = r1-r0;
            msd_liquid(stepSize, j) = msd_liquid(stepSize, j) + dot(d, d); % Add displacement of the particle
        end
    end
end


close(h);
close(h2);
disp('Beep!');
beep;

plot_msd;