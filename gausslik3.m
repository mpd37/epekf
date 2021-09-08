function [logZ dlogZdm dlogZdS dlogZdx] = gausslik3(m, S, x)

% inputs:
% m
% S
% x
%
% returns logZ, dlogZdm, dlogZdS, dlogZdx
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth 

D = size(S,1); % dimensionality of the distribution

diff_vec = bsxfun(@minus,x,m);
x2 = S\diff_vec;

z = (2*pi)^(D/2)*sqrt(det(S)); % inverse normalization constant


logZ = -sum(diff_vec.*x2,1)/2 - log(z); % log likelihood

if nargout > 1
  dlogZdm = (x-m)'/S; % OK
  dlogZdS = dlogZdm'*dlogZdm/2 - pinv(S)/2; % OK
  dlogZdx = -dlogZdm; % OK
end
