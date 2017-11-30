Zu = 2 % Unscreened
Zvov = 27/16 %variationally optimized value
p = @(r, Z) Z.^3*4*r.^2 .*exp(-2.*Z.*r);

r = linspace(0, 4);

figure;
hold on
plot(r, p(r, Zu));
plot(r, p(r, Zvov));

legend({'$p(r), Z = 2$', '$p(r), Z = \frac{27}{16}$'}, 'Interpreter', 'latex');
title('P(r)')