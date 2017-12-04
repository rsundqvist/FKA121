function localEnergy = computeLocalEnergy(r_vec1,r_vec2,alpha)

r1 = norm(r_vec1);
r2 = norm(r_vec2);
r12 = norm(r_vec1 - r_vec2);

r1Hat = r_vec1/r1;
r2Hat = r_vec2/r2;
dotTerm = dot(r1Hat-r2Hat,r_vec1 - r_vec2)/(r12*(1 + alpha*r12)^2);

localEnergy = -4 + dotTerm - 1/(r12*(1 + alpha*r12)^3) - 1/(4*(1 + alpha*r12)^4) + 1/r12;

end