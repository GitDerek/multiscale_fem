function [phi,dphids,dphidt] = Q1basis(s,t)

% definition of Q1 shape function
phi(1) = 0.25*(s-1.)*(t-1.);
phi(2) = -0.25*(s+1.)*(t-1.);
phi(3) = 0.25*(s+1.)*(t+1.);
phi(4) = -0.25*(s-1.)*(t+1.);
%
dphids(1) = 0.25*(t-1.);
dphids(2) = -0.25*(t-1.);
dphids(3) = 0.25*(t+1.);
dphids(4) = -0.25*(t+1.);
%
dphidt(1) = 0.25*(s-1.);
dphidt(2) = -0.25*(s+1.);
dphidt(3) = 0.25*(s+1.);
dphidt(4) = -0.25*(s-1.);
end