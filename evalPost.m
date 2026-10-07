function postker = evalPost(mu,para,data,i)
% Function EVALPOST
%
% Purpose:    Evaluate log posterior kernel
%
% Format:     postker = evalPost(mu,para,data,i)
%
% Input:      mu          current state mean
%             para        parameter structure
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
%             i           current state
%
% Output:     postker     (-)log posterior kernel
%
% Written by Fei Tan, Saint Louis University
% Updated: September 1, 2022

%% -------------------------------------------
%          Evaluate Posterior Kernel
%---------------------------------------------

% Evaluate prior
[ybar,logprior] = ssEqm(mu);

% Evaluate posterior
if isinf(logprior)
    postker = 1e10;
else
    n = size(data.Y,2);
    nlag = size(data.X,2)/n;
    ind = find(para.S==i);   % non-empty
    T = length(ind);
    Y = data.Y(ind,:)-repmat(ybar,T,1);
    X = data.X(ind,:)-repmat(ybar,T,nlag);
    U = Y-X*para.Phi(:,:,i);
    loglik = -n*T/2*log(2*pi)-T/2*logdet(para.Sig(:,:,i),'chol')-trace(para.Sig(:,:,i)\U'*U)/2;
    postker = -(loglik+logprior);
end

%-------------------- END --------------------
