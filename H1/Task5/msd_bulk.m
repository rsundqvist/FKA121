%% Plot vol
close all;

%loadAllpos;

dt = 1;
numRows = size(data, 1);

beep on;

h = waitbar(0, 'Progress - overall');
h2 = waitbar(0, 'Progress - particle');
maxOffset = 1; % Step size if amplified by ir in MD_main.c!

ir = 5;
N = 256;

solid_i0 = 1;
msd_solid = zeros(maxOffset, numRows);
for offset=1:maxOffset
    waitbar(0.5*offset/maxOffset, h);
    for i=1:N
        waitbar(i/N, h2);
        for t = 1:numRows-offset
            i1 = 1+(i-1)*3; i2 = 3+(i-1)*3; % x to z for the particle
            r1 = pos(t+offset,i1:i2); % Position after stepSize iterations
            r0 = pos(t,i1:i2); % Base position
            r0 = pos(solid_i0, i1:i2);
            d = r1-r0;
            msd_solid(offset, t) = msd_solid(offset, t) + dot(d, d); % Add displacement of the particle
        end
    end
end

liquid_i0 = 5000;
msd_liquid = zeros(maxOffset, numRows);
for offset=1:maxOffset
    waitbar(0.5 + offset/maxOffset, h);
    for i=1:N
        waitbar(i/N, h2);
        for t = 1:numRows-offset
            i1 = 1+(i-1)*3; i2 = 3+(i-1)*3; % x to z for the particle
            r1 = pos2(t+offset,i1:i2); % Position after stepSize iterations
            r0 = pos2(t,i1:i2); % Base position
            r0 = pos2(liquid_i0,i1:i2);
            d = r1-r0;
            msd_liquid(offset, t) = msd_liquid(offset, t) + dot(d, d); % Add displacement of the particle
        end
    end
end



close(h);
close(h2);
disp('Beep!');
beep;

plot_msd;