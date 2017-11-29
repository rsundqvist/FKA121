
figure
hold on;

series = msd_solid(stepSize,:)/N;
%series = series.*series;
plot(series, 'DisplayName', ['Solid: stepSize = ' num2str(stepSize*ir)]);
    
hold off;
legend(gca, 'show');
title('Title goes here')