clear all;
clc;
% Compute the norm for our trial wave function
alpha = 0.1;
psiT = @(x1,y1,z1,x2,y2,z2) exp(-4*(x1+x2+y1+y2+z1+z2))* ...
    exp(sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)/1 + 0.1*sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2));

start = [0,0,0,0,0,0];
nsamples = 10^4;

% bullshit
proppdf = @(x,y)gampdf(x,floor(alpha),floor(alpha)/alpha);
proprnd = @(x)sum(...
              exprnd(floor(alpha)/alpha,floor(alpha),1));

% Metropolis
[smpl,accept] = mhsample(start,nsamples,'pdf',psiT,'proppdf',proppdf, 'proprnd',proprnd)