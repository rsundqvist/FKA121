N  = 10000000;  
data = zeros(1,N);
for i=1:N
    data(i) = rand;
end

h = histogram(data,'normalization','pdf');

b = h.BinWidth;
a = h.BinCounts;
c = h.BinEdges;
sum(a*b)

figure(2)
plot(c(2:end),a/sum(a*b))