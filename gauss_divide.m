% function [m v] = gauss_divide(m1, v1, m2, v2)
%
% v = 1/(1/v1 - 1/v2);
% m = v*(m1/v1 - m2/v2);
%
% Marc Deisenroth, 2012-05-28


function [m S dmda dmdA dSda dSdA] = gauss_divide(a, A, b, B)

D = size(A,1);

if any(any(B==Inf)) || any(any(isnan(B))) % some funny cases
  S = A;
  m = a;
else
  ABB = pinv(A - B)*B;
  S = -A*ABB; % new covariance matrix: (A\inv - B\inv)\inv
  S = (S+S')/2;
  m = S*(pinv(A)*a - pinv(B)*b); % S*(A\inv a - B\inv b)
end

if nargout > 2
  dSda = zeros(D,D,D);
  dSdA = zeros(D,D,D,D);
  
  dmda = S*pinv(A); % OK
  dmdb = -S*pinv(B);
  
  iAB = pinv(A-B);
  iA = pinv(A);
  for i = 1:D
    for j = 1:D
      
      Eij = zeros(D); Eij(i,j) = 1;
      dSdA(:,:,i,j) = -Eij*ABB + (A*iAB(:,i))*(iAB(j,:)*B); % OK
      dmdA2(:,i,j) = dSdA(:,:,i,j)*(pinv(A)*a - pinv(B)*b);
      dmdA_inner(:,i,j) =  (Eij*a);
    end
  end
  
  dmdA1 = -etprod('134',etprod('134',S*iA,'24',dmdA_inner,'123'),'123',iA,'24'); % OK (derivative of the bracket wrt A);
  
  dmdA = dmdA1 + dmdA2;
  
end

