function msd = computeMSD(pos,tmax,N)

maxStep = 10000
msd = zeros(1,maxStep);


%for step=1:maxStep
%    for t=1:tmax-step
%        tmpSumP = 0; % sum over all particle displacements at given t
%        for i=1:N
%             
%             i1 = 1+(i-1)*3;
%             i2 = 3+(i-1)*3;
%             
%             r1 = pos(t,i1:i2);
%             r2 = pos(t+step,i1:i2);
%             tmpSumP = tmpSumP + norm(r2-r1)^2;
%         end
%         msd(step) = msd(step) + tmpSumP/N;
%     end
%     msd(step) = msd(step)/(tmax-step);
%     disp(step);
% end

for step=1:maxStep
    tmpSum = 0;
    for i=1:N
             i1 = 1+(i-1)*3;
             i2 = 3+(i-1)*3;
             
             r1 = pos(1,i1:i2);
             r2 = pos(1+step,i1:i2);
             tmpSum = tmpSum + norm(r1-r2)^2;
    end
    msd(step) = tmpSum/N;
end

end