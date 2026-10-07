function [predS,updtS,loglik] = evalLik(para,data)
% Function EVALLIK
%
% Purpose:    Evaluate VAR log likelihood
%
% Format:     [predS,updtS,loglik] = evalLik(para,data)
%
% Input:      para        parameter structure
%
%                         .mu(:,i)    - state-i mean, i = 1,...,ns
%                         .Phi(:,:,i) - state-i coefficient, i = 1,...,ns
%                         .Sig(:,:,i) - state-i covariance, i = 1,...,ns
%                         .tp(i)      - state-i transition prob, i = 1,...,ns
%                         .S(t)       - state at time t, t = 1,...,T
%
%             data        data structure
%
%                         .Y(t,:)     - current observation at time t, t = 1,...,T
%                         .X(t,:)     - lagged observations at time t, t = 1,...,T
%
% Output:     predS(:,t)  P[s(t)|Y(:,1:t-1)], t = 1,...,T
%             updtS(:,t)  P[s(t)|Y(:,1:t)], t = 1,...,T
%             loglik(t)   log p(y(t)|Y(:,1:t-1)), t = 1,...,T
%
% Written by Fei Tan, Saint Louis University
% Updated: September 1, 2022

%% -------------------------------------------
%            Likelihood Evaluation
%---------------------------------------------

% Initialization
[T,n] = size(data.Y);        % numbers of periods & variables
nlag = size(data.X,2)/n;     % number of VAR lags
ns = length(para.tp);        % number of states
predS = zeros(ns,T+1);       % predictive state prob
predS(1,1) = 1;
updtS = zeros(ns,T);         % updated state prob
loglik = zeros(T,1);         % period log likelihood

% Evaluate likelihood
for t = 1:T
    % Period-t updated state prob
    for i = 1:ns
        [ybar,~] = ssEqm(para.mu(:,i));
        yt = data.Y(t,:)-ybar;
        xt = data.X(t,:)-repmat(ybar,1,nlag);
        updtS(i,t) = exp(mvt_pdf(yt,xt*para.Phi(:,:,i),para.Sig(:,:,i),inf))*predS(i,t);
    end
    loglik(t) = sum(updtS(:,t));
    updtS(:,t) = updtS(:,t)/loglik(t);
    
    % Period-(t+1) predictive state prob
    for i = 1:ns
        if i==1
            predS(i,t+1) = para.tp(i)*updtS(i,t);
        else
            predS(i,t+1) = (1-para.tp(i-1))*updtS(i-1,t)+para.tp(i)*updtS(i,t);
        end
    end
end
predS = predS(:,1:end-1);
loglik = log(loglik);

%-------------------- END --------------------
