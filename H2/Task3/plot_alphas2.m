alphas = 0.05:0.005:0.25;
sz = length(alphas);

avg = zeros(sz,1);
err = zeros(sz,1);
N = 5000;

for i = 1:sz
    name = ['alphas' num2str(i-1) '.dat'];
    data = load(name);
    
    avg(i) = mean(data(:,1));
    err(i) = mean(data(:,2))/sqrt(N);
end

errorbar(alphas, avg, err);