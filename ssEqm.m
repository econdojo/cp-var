function [ybar,logprior] = ssEqm(mu)
% Function SSEQM
%
% Purpose:    Impose steady-state equilibrium
%
% Format:     [ybar,logprior] = ssEqm(mu)
%
% Input:      mu        [gbar pbar sbar dRbar]
%
% Output:     ybar      [gbar pbar Rbar sbar dRbar]
%             logprior  log prior
%
% Written by Fei Tan, Saint Louis University
% Updated: September 1, 2022

%% -------------------------------------------
%              S.S. Equilibrium
%---------------------------------------------

% s.s. eqm
gbar = mu(1);
pbar = mu(2);
sbar = mu(3);
dRbar = mu(4);
Rbar = pbar+gbar+sbar;
ybar = [gbar pbar Rbar sbar dRbar];

% log prior
if pbar<=0 || Rbar<=0 || dRbar<=0
    logprior = -inf;
else
    logprior = prior_pdf(gbar,0.4,0.2,'t')...
          +prior_pdf(pbar,0.8,0.5,'G')...
          +prior_pdf(sbar,0,2,'t')...
          +prior_pdf(dRbar,0.5,0.2,'G');
end

%-------------------- END --------------------