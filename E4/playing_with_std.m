
x = zeros(10,100);

for i=2:1000
   for j=1:10
      x(j,i) = x(j,i-1) + rand; 
   end
end

g = zeros(1,1000);

for i=1:1000
   g(i) = std(x(:,i));
end

figure(555)
plot(g)