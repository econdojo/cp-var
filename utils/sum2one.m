function v = sum2one(v)
% Function SUM2ONE
%
% Purpose:    Normalize a vector to sum to one
%
% Format:     v = sum2one(v)
%
% Input:      v         vector of non-negative numbers
%
% Output:     v         normalized vector
%
% Written by Fei Tan, Saint Louis University
% Updated: September 1, 2022

% Normalization
realsmall = 1e-6;
v = v/sum(v);
if sum(v)~=1
    [~,i] = max(v);
    v(i) = v(i)-realsmall;
    v(end) = 1-sum(v(1:end-1));
end

%-------------------- END --------------------