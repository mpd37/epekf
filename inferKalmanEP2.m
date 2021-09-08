function out = inferKalmanEP2(x, y, m0, S0, FF, GG, Q, R, u)

% Function that implements the extended Kalman smoother as an iterative message
% passing algorithm by means of EP
%
%
% Input arguments:
%
% x     sequence of latent states (if available)
% y     sequence of observations
% m0    prior mean over latent state
% S0    prior covariance over latent state
% FF    function handle to transition model
% GG    function handle to measurement model
% Q     covariance matrix of system noise
% R     covariance matrix of measurement noise
% u     sequence of control signals (if applicable)
%
% Copyright (C) 2012-2013 by 
% Marc Deisenroth and Shakir Mohamed
%
% Last modified: 2013-07-18

EPiter = 100;
display = true;

%% Initialization

% Check if there is a control signal
control = false; dimU = 0;
if ~isempty(u)
  control = true;
  dimU = size(u,1);
end;

d = size(S0,1); % (latent) state dimension
e = size(y, 1); % dimension of observation
n = size(y, 2); % length of time series

mx = cell(1,n+1);
vx = cell(1,n+1);

% initialize the messages and the marginals
for i = 1:n+1
  fw(i).mean = zeros(d,1); fw(i).cov = Inf*eye(d); fw(i).s = 1;
  bw(i).mean = zeros(d,1); bw(i).cov = Inf*eye(d); bw(i).s = 1;
  meas(i).mean = zeros(d,1); meas(i).cov = Inf*eye(d); meas(i).s = 1;
  mx{i} = zeros(d,1); vx{i} = 100000*eye(d);
end;

mx_old = []; vx_old = [];

% prior distribution put into the first variable node
mx{1} = m0;  vx{1} = S0;

%% Run EP
for iter = 1:EPiter
  
  % Forward Sweep
  for i = 1:n % go step-by-step through the time series
    fprintf('\rEPiter #%i/%i   Node %i/%i', iter, EPiter, i, n);
    %% Foward messages
    if i > 1
      % Cavity distribution of curr node (distr excl. fwd msg)
      [mu_xfwd, Sigma_xfwd] = gauss_divide(mx{i}, vx{i}, fw(i).mean, fw(i).cov);
      
      % Cavity distribution of prev node (distr excl. back msg)
      [mu_xback, Sigma_xback] = gauss_divide(mx{i-1}, vx{i-1}, bw(i-1).mean, bw(i-1).cov);
      
      
      % Update the forward messages
      % linearize
      if control
        [Mx, A] = jaccsd(FF, mu_xback, u(:,i-1));
      else
        [Mx, A] = jaccsd(FF, mu_xback);
      end
      
      fw(i).mean = Mx; fw(i).cov = Q + A*Sigma_xback*A';
      [logZi, dlogZidMx, dlogZidSx] = gausslik3(fw(i).mean, fw(i).cov + Sigma_xfwd, mu_xfwd);
      
      % Update marginal using updated message
      if iter == 1 % Sigma_xfwd = Inf
        [mx{i}, b] = gm_w_prec(mu_xfwd, zeros(d), fw(i).mean, pinv(fw(i).cov));
        vx{i} = pinv(b);
      else
        [mx{i}, vx{i}] = gm(mu_xfwd, Sigma_xfwd, fw(i).mean, fw(i).cov);
      end
      
      if det(vx{i}) < 0
        fprintf('Forward message: negative variance in posterior marginal ...\n');
        fprintf('Check wether we ended up outside the training data!');
      end
    end;
    
    %% ---- Measurement messages
    % compute cavity distribution
    [mu_xup, Sigma_xup ] = gauss_divide(mx{i}, vx{i}, meas(i).mean, meas(i).cov);
    
    if any(eig(Sigma_xup)<0);
      fprintf('break the loop (measurement message)\n'); continue;
    end
    
    % Compute predictive distr using GP with uncertain inputs
    [mz, C] = jaccsd(GG, mu_xup);
    sz = C*Sigma_xup*C' + R;
    dMzdmu_xup = C;
    
    % If we have an injective mapping, go via the derivatives
    if d <= e
      [logZi dlogZidMz dlogZidSz] = gausslik3(mz, sz, y(:,i));
      dlogZidmu_xup = dlogZidMz*dMzdmu_xup;
      
      dSzdSigma_xup = zeros(e,e,d,d);
      for k = 1:d
        for l = 1:d
          dSzdSigma_xup(:,:,k,l) = dMzdmu_xup(:,k)*dMzdmu_xup(:,l)'; % assume constant C
        end
      end
      dlogZidSigma_xup = reshape(dlogZidSz(:)'*reshape(dSzdSigma_xup, e*e, d*d), d,d);
      
      
      % Update Marginal
      [mx{i}, vx{i}] = moment_matching_gaussian(mu_xup, Sigma_xup, dlogZidmu_xup, dlogZidSigma_xup);
      
      % Update the message using derivatives of logZi
      [meas(i).mean, meas(i).cov, meas(i).s] ...
        = site_update_gaussian(mu_xup, Sigma_xup, exp(logZi), dlogZidmu_xup, dlogZidSigma_xup);
      
    else % follow Yuan Qi's description (no derivatives needed)
      C = dMzdmu_xup';
      Vxz = Sigma_xup*C;
      K = Vxz/sz;
      mx{i} = mu_xup + K*(y(:,i) - mz);
      vx{i} = Sigma_xup - K*Vxz';
      [meas(i).mean, meas(i).cov] = gauss_divide(mx{i}, vx{i}, mu_xup, Sigma_xup);
    end
  end;
  
  
  
  
  %% --- Backward messages
  for i = n-1:-1:1
    fprintf('\rEPiter #%i   Node %i/%i', iter, i, n);
    
    % Cavity Distributions
    [mu_xback, Sigma_xback] = gauss_divide(mx{i}, vx{i}, bw(i).mean, bw(i).cov);
    [mu_xfwd, Sigma_xfwd] = gauss_divide(mx{i+1}, vx{i+1}, fw(i+1).mean, fw(i+1).cov);
    
    if any(eig(Sigma_xback)<0) | any(eig(Sigma_xfwd)<0)
      fprintf('break the loop (backward message)\n'); continue;
    end;
    
    % Compute predictive distribution
    if control
      [Mx, A] = jaccsd(FF, mu_xback, u(:,i));
    else
      [Mx, A] = jaccsd(FF, mu_xback); % "exact" prediction
    end
    Sx = A*Sigma_xback*A' + Q;
    dMxdmxni = A;
    
    
    if true
      [logZi dlogZidMx dlogZidSx] = gausslik3(Mx, Sx + Sigma_xfwd, mu_xfwd);
      
      % chain-rule
      %       dlogZidMxni = etprod('2',dlogZidMx(:), '1', dMxdmxni(:,1:d), '12')';
      dlogZidMxni = dlogZidMx*dMxdmxni(:,1:d);
      dSxdsxni = zeros(d,d,d,d);
      for k = 1:d
        for l = 1:d
          dSxdsxni(:,:,k,l) = dMxdmxni(:,k)*dMxdmxni(:,l)';
        end
      end
      %       dlogZidVxni = etprod('34',dlogZidSx,'12',dSxdsxni,'1234');
      dlogZidVxni = reshape(dlogZidSx(:)'*reshape(dSxdsxni, d*d, d*d),d,d);
      
      % Update Marginal
      [mx{i}, vx{i}] = moment_matching_gaussian(mu_xback(1:d), Sigma_xback(1:d,1:d), dlogZidMxni, dlogZidVxni);
       % Update the backward message
      [bw(i).mean, bw(i).cov, bw(i).s] ...
      = site_update_gaussian(mu_xback(1:d), Sigma_xback(1:d,1:d), exp(logZi), dlogZidMxni, dlogZidVxni);
    else
      % Qi's approach (derivative-free)
      Vxx = Sigma_xback(1:d,1:d)*dMxdmxni(:,1:d)'; % cross-covariance
      J = Vxx(1:d,:)/Sx;
      % correct for error between predicted mean (from cavity distrbution and marginal mean)
      mx{i} = mu_xback(1:d,:) + J(1:d,:)*(mx{i+1} - Mx);
      vx{i} = Sigma_xback(1:d,1:d) + J(1:d,:)*(vx{i+1} - Sx)*J';
      
      [bw(i).mean, bw(i).cov] = gauss_divide(mx{i}, vx{i}, mu_xback, Sigma_xback);
    end

  end
  
  % Collect results
  mm{iter} = cell2mat(mx);
  vv{iter} = vx;
  
  
  if exist('x', 'var')
    for nn = 1:n
      nll(nn) = -gausslik3(mm{iter}(:,nn), vv{iter}{nn}, x(:,nn));
    end
  fprintf('   NLL (x), %g\n', mean(nll));  
  end
  
  
  
  % check for convergence
  if iter > 1
    meanCriterion = norm(cell2mat(mx)-cell2mat(mx_old))/n;
    varCriterion = norm(cell2mat(vx)-cell2mat(vx_old))/n;
    if (meanCriterion < 1e-6) && (varCriterion < 1e-6)
      fprintf('EP converged after %d iterations \n', iter);
      break;
    end;
  end;
  mx_old = mx; vx_old = vx;
  
end;

% Do a check
tmp = cellfun(@(x)sum(diag(x)<0),vv{end});
if sum(tmp), fprintf('\nWarning! Negative variances in final marginals \n'); end;

% Return results
out.mean = mm;
out.cov = vv;
out.messages.fw = fw;
out.messages.bw = bw;
out.messages.meas = meas;

