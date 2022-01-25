% Demo script for the EP version of the extended Kalman filter
%
% Example: pendulum motion
%
% Copyright (C) 2012-2022 by 
% Marc Deisenroth and Shakir Mohamed
%

clear all; close all;

randn('state',1); rand('state', 2);

n = 50; % length of time series

generate_data_pendulum;



%% Inference

out = inferKalmanEP2(x, y, m0, S0, FF, GG, Q, R, u);


%% Compute Error Metrics
[~,T] = size(y);
xAvail = ~isempty(x);
nlpX = NaN; maeX = NaN;


% compute the negative log-likelihood
for i = 1:size(out.mean,2) % for all iters
  for t = 1:T % for each observations
    if xAvail
      nllpX_tmp(t) = -gausslik3(out.mean{i}(:,t), out.cov{i}{t}, x(:,t));
    end;
  end;
  
  if xAvail
    nlpX(i) = mean(nllpX_tmp);
  end;
end;


figure;
hold on
plot(1:length(nlpX), nlpX, '-o')
xlabel('EP iteration');
ylabel('Neg. log likelihood in x')
axis tight