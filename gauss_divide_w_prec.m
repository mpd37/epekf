function [m iS] = gauss_divide_w_prec(a, iA, b, iB)

% function [m iS] = gauss_divide_w_prec(a, iA, b, iB)
% 
% divides first Gaussian by second Gaussian
% 
% takes means and two Precisions
% returns new mean and new precision; OK
%
% Marc Deisenroth, 2012-05-21

iS = iA - iB;

m = pinv(iS)*(iA*a - iB*b); % S*(A\a - B\b)


