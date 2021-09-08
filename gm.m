function [m, S, z] = gm(a, A, b, B)

% multiplies two Gaussians together, whose means are given by a, b, and
% whose covariances are given by A, B.
%
% with 4 inputs: call gm(a, A, b, B)
%
% if only 3 inputs are given (a,A,b), it is assumed that B=A;
%
% if only 2 inputs are given (a,b), it is assumed that A=I=B
%
% returns:
% m: mean of resulting Gaussian
% S: covariance of resulting Gaussian
% z: normalization constant of resulting Gaussian
%
% (C) 2008-2013 by Marc Deisenroth
% Last modified: 2010-07-30


a = a(:); % enforce column vectors
D = length(a); % dimension

switch nargin
  
  case 4 % multiply two different Gaussians
    AB = A+B;
    %     L = chol(AB)';              % cholesky factorization of the covariance
    
    S = A*((AB)\B); % covariance matrix
    m = B*((AB)\a) + A*((AB)\b); % mean
    
    
    
    if nargout == 3
      alpha = AB\(a-b);             % precomputation
      z = exp(-0.5*sum((a-b).*alpha,1))./((2*pi)^(0.5*D)*sqrt(det(AB)));
      keyboard
    end
    
  case 3 % assume 2 Gaussians with the same covariance (A = B)
    b = b(:);
    AB = 2*A;
    L = chol(AB)';
    alpha = L\(a-b);
    
    S = A/2;
    m = (a+b)/2;
    
    if nargout == 3
      z = exp(-0.5*sum(alpha.^2,1))./((2*pi)^(0.5*D)*prod(diag(L)));
    end
    
  case 2 % assume 2 unit-covariance normals with different means
    b = A(:);
    alpha = (a-b)./sqrt(2);
    
    S = eye(D)/2;
    m = (a+b)/2;
    
    if nargout == 3
      z = exp(-0.5*sum(alpha.^2,1))./((4*pi)^(0.5*D));
    end
    
end

S = (S+S')/2;