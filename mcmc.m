% MAIN function: implement MCMC algorithm
% Written by Fei Tan, Saint Louis University
% Updated: September 1, 2022

clear
close all
clc
addpath(genpath('utils'));

%% -------------------------------------------
%               MCMC Algorithm
%---------------------------------------------

% Initialization
M = 11000;
burn = 1000;
nlag = 4;
np = 4;
ns = 3;
tp_a = 6;
tp_b = 1;
pdof = 15;
Y = importdata('data.txt');
[T1,n] = size(Y);
T1 = T1-nlag;
chain.mu = zeros(np,ns,M);
chain.Phi = zeros(n*nlag,n,ns,M);
chain.Sig = zeros(n,n,ns,M);
chain.tp = zeros(M,ns-1);
chain.S = zeros(M,T1);
chain.predS = zeros(ns,T1,M);
para.mu = zeros(np,ns);
para.Phi = zeros(n*nlag,n,ns);
para.Sig = zeros(n,n,ns);
para.tp = zeros(1,ns);
para.S = zeros(1,T1);
ind = cell(ns,1);
T2 = zeros(ns,1);
rej1 = zeros(ns,1);
rej2 = 0;

% Construct VAR data
data.Y = Y((nlag+1):end,:);  % exclude initial lags
data.X = zeros(T1,n*nlag);   % regressors
for k = 1:nlag
    data.X(:,((k-1)*n+1):(k*n)) = Y((nlag-k+1):(nlag-k+T1),:);
end

% Current parameter draw
para.tp(1:ns-1) = tp_a/(tp_a+tp_b);
para.tp(end) = 1;
para.S(1:72) = 1;       % 1963 - 1980
para.S(73:152) = 2;     % 1981 - 2000
para.S(153:end) = ns;   % 2001 - 2021
for i = 1:ns
    ind{i} = find(para.S==i);
    T2(i) = length(ind{i});
    ybar = mean(data.Y(ind{i},:));
    para.mu(:,i) = ybar([1 2 4 5])';
    [ybar,~] = ssEqm(para.mu(:,i));
    Y = data.Y(ind{i},:)-repmat(ybar,T2(i),1);
    X = data.X(ind{i},:)-repmat(ybar,T2(i),nlag);
    para.Phi(:,:,i) = (X'*X)\(X'*Y);
    para.Sig(:,:,i) = (Y'*Y-(Y'*X)*para.Phi(:,:,i))/T2(i);
    R = cholmod(para.Sig(:,:,i)); para.Sig(:,:,i) = R'*R;
end

% MCMC
tic; progressbar('MCMC in Progress')
for iter = 1:M
    % Sample parameters, conditional on states
    for i = 1:ns
        % mu
        [~,xopt,~,Hopt,~,~,~] = csminwel('evalPost',para.mu(:,i)',eye(np)*1e-3,[],1e-5,100,1,para,data,i);
        mu = mvt_rnd(xopt,Hopt,pdof,1);
        pkNext = -evalPost(mu,para,data,i);
        pkLast = -evalPost(para.mu(:,i)',para,data,i);
        r = exp(pkNext-pkLast+mvt_pdf(para.mu(:,i)',xopt,Hopt,pdof)-mvt_pdf(mu,xopt,Hopt,pdof));
        if rand > min([r 1])
            chain.mu(:,i,iter) = para.mu(:,i);
            rej1(i) = rej1(i)+1;
        else
            chain.mu(:,i,iter) = mu';
            para.mu(:,i) = mu';
        end
        
        % Phi, Sigma
        [ybar,~] = ssEqm(para.mu(:,i));
        Y = data.Y(ind{i},:)-repmat(ybar,T2(i),1);
        X = data.X(ind{i},:)-repmat(ybar,T2(i),nlag);
        Phi = (X'*X)\(X'*Y);
        Sig = (Y'*Y-(Y'*X)*Phi)/T2(i);
        R = cholmod(Sig); Sig = R'*R;
        chain.Sig(:,:,i,iter) = iwishrnd(T2(i)*Sig,T2(i)+n-1);
        chain.Phi(:,:,i,iter) = reshape(mvt_rnd(reshape(Phi,n*nlag*n,1)',kron(chain.Sig(:,:,i,iter),inv(X'*X)),Inf,1)',n*nlag,n);
        para.Sig(:,:,i) = chain.Sig(:,:,i,iter);
        para.Phi(:,:,i) = chain.Phi(:,:,i,iter);
        
        % transition probabilities
        if i<ns
            chain.tp(iter,i) = betarnd(tp_a+T2(i)-1,tp_b+1);
            para.tp(i) = chain.tp(iter,i);
        end
    end
    
    % Sample states
    [chain.predS(:,:,iter),updtS,~] = evalLik(para,data);
    S = para.S;
    for t = T1-1:-1:2
        if para.S(t+1)>1
            p1 = updtS(para.S(t+1)-1,t)*(1-para.tp(para.S(t+1)-1));
            p2 = updtS(para.S(t+1),t)*para.tp(para.S(t+1));
            if rand <= p1/(p1+p2)
                para.S(t) = para.S(t+1)-1;
            else
                para.S(t) = para.S(t+1);
            end
        else
            para.S(t) = 1;
        end
    end
    
    ind2 = ind; T3 = T2;
    for i = 1:ns
        ind{i} = find(para.S==i);
        T2(i) = length(ind{i});
    end
    if min(T2)<10
        para.S = S;
        ind = ind2;
        T2 = T3;
        rej2 = rej2+1;
    end
    chain.S(iter,:) = para.S;
    
    progressbar(iter/M)
    clc
end

%% -------------------------------------------
%             Posterior Analysis
%---------------------------------------------

rej1 = rej1/M;
M = M-burn;
chain = chain((burn+1):end,:);
stat = PostStat(chain,para(:,1));
num = num((burn+1):end);
den = den((burn+1):end);
logmlik = pk_opt-log(mean(num)/mean(den));
H = [num den];
H = H-repmat(mean(H),M,1);
dev = [1./mean(num) -1./mean(den)];
nse = sqrt(dev*NeweyWest(H)*dev'/M);

time = datestr(toc/(24*60*60),'HH:MM:SS');
fprintf('\n');
fprintf('Number of draws = %d after %d burn-in\n',M,burn);
fprintf('Elapsed time = %s [hh:mm:ss]\n',time);
fprintf('Rejection rate  =  %.1f %%\n',rej1*100);
fprintf('Log marginal likelihood  =  %.3f (n.s.e. %.4f)\n\n',logmlik,nse);
save([mod filesep 'mcmc.mat'],'chain','stat','num','den','logmlik','nse','xopt','Hopt');

if np<=4
    dim1 = 1;
elseif np>=5 && np<=8
    dim1 = 2;
elseif np>=9 && np<=12
    dim1 = 3;
elseif np>=13
    dim1 = 4;
end
dim2 = ceil(np/dim1);
for k = 1:np
    subplot(dim1,dim2,k);
    pd = fitdist(chain(:,k),'Kernel','Kernel','epanechnikov');
    x = linspace(stat(k,2)-0.2,stat(k,3)+0.2);
    d = pdf(pd,x);
    plot(x,d,'--r','LineWidth',1)
    xlabel(para{k,1})
end
saveas(gcf,[mod filesep 'Fig_postpdf.fig']);
