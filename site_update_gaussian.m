function [mi, vi, si] = site_update_gaussian(mni, vni, Zi, dlogZdm, dlogZdv)
% Update a site distribution based on the cavity distribution, the
% log-partitition function, and its derivatives
%
%
% Copyright (C) 2012-2013 by 
% Marc Deisenroth and Shakir Mohamed
%
% Last modified: 2013-07-18

dlogZdm = dlogZdm(:)'; % make sure it's a row vector

d = size(vni,2);

vi = pinv(dlogZdm'*dlogZdm - 2*dlogZdv) - vni; % OK
mi = mni + (dlogZdm'*dlogZdm - 2*dlogZdv)\dlogZdm'; % OK
si = Zi*sqrt(det(eye(d) + vni/vi))*exp(dlogZdm*((dlogZdm'*dlogZdm - 2*dlogZdv)\dlogZdm')/2); % OK
si = max(si, 1e-8);

