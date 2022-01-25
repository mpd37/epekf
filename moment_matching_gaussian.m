function [mx, vx] = moment_matching_gaussian(mxni, vxni, dlogZidmxni, dlogZidvxni)


% EP updates of the moments of the factored distribution
%
% from T. Minka: "EP: A quick reference", 2008
%
% Marc Deisenroth, 2012-02-13

dlogZidmxni = dlogZidmxni(:)'; % make sure it's a row vector

mx = mxni + vxni*dlogZidmxni'; % soft update
vx = vxni + vxni*(2*dlogZidvxni - dlogZidmxni'*dlogZidmxni)*vxni;