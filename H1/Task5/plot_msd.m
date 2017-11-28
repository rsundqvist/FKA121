
figure
hold on;
for stepSize=1:maxStepSize
    series = msd_solid(stepSize,:)/N;
    %series = series.*series;
    plot(series, 'DisplayName', ['Solid: stepSize = ' num2str(stepSize*ir)]);
    
    series = msd_liquid(stepSize,:)/N;
    %series = series.*series;
    plot(series, 'DisplayName', ['Liquid: stepSize = ' num2str(stepSize*ir)]);
end
hold off;
legend(gca, 'show');
title('Title goes here')