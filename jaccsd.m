function [z,A]=jaccsd(fun,x,u)
% JACCSD Jacobian through finite differences
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
% Copyright (C) 2012-2013 by 
% Marc Deisenroth and Shakir Mohamed
%
% Last modified: 2013-07-18

if exist('u','var')
  z = fun(x,u);
else
  z = fun(x);
end
n=numel(x);
m=numel(z);
A=zeros(m,n);
h = 1e-4;
for k=1:n
    x1=x;
    x1(k)=x1(k) + h;
    if exist('u','var')
      A(:,k)=(fun(x1,u)-z)/h;
    else
      A(:,k)=(fun(x1)-z)/h;
    end
end