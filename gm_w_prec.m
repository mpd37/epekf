function [m iS] = gm_w_prec(a, iA, b, iB)

% 
% multiplies first Gaussian by second Gaussian
% 
% takes means and two Precisions
% returns new mean and new precision; OK
%
% Marc Deisenroth, 2010-11-17

iS = iA + iB;
m = iS\(iA*a + iB*b); % S*(A\inv a - B\inv b)
