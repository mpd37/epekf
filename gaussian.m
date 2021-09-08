function x = gaussian(m, S, n)

% Generate samples from a Gaussian
%
% Copyright (C) 2007-2013 by 
% Marc Deisenroth and Shakir Mohamed
%
% Last modified: 2013-07-18

if nargin < 3, n = 1; end

x = bsxfun(@plus, m(:), chol(S)'*randn(length(m),n));
