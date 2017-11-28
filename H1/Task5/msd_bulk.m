%% Plot vol
close all;

%loadAllpos;
%t = data(:,1);

dt = 1;
len = length(data);

beep on;

h = waitbar(0, 'Progress');
diff = zeros(len, 256);
diff2 = zeros(len, 256);
%%
maxStepSize = 5
timeSteps = len;
msd = zeros(1,maxStepSize);
for stepSize=1:maxStepSize
    for i=1:256
        for j = 1: timeSteps - stepSize
            i1 = 1+(r-1)*3;
            i2 = 3+(r-1)*3;
            r2 = pos(j+stepSize,i1:i2);
            r2 = pos(j,i1:i2);
        end
    end
end

%disp('Begin work');
%for i = 1:len-dt
%    waitbar(i/len, h)
%    for j = 1:3:256
%        X0 = pos(i,j:j+2);
%        X = pos(dt+i,j:j+2);
%        
%        
%    end
%    
%    diff(i) = norm(pos(i+dt,:) - pos(i,:))/256;
%    diff2(i) = norm(pos2(i+dt,:) - pos2(i,:))/256;
%end






close(h);
disp('Beep!');
beep;

hold on;
plot(t, diff);
plot(t, diff2);

xlabel('Time [ps]')
ylabel('msd')
legend('Liquid', 'Solid');

hold off;
