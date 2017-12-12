x = linspace(0,10);
y = sin(3*x);
yyaxis left
plot(x,y)

z = sin(3*x).*exp(0.5*x);
yyaxis right
plot(x,z)
ylim([-150 150])