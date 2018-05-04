function [mESS,Sigma,b] = multiESS(X,Sigma,b,Noffsets,Nb)
%MULTIESS Compute multivariate effective sample size of Markov chain.
%   MESS = MULTIESS(X) computes effective sample size MESS of single Markov 
%   chain X, using the multivariate dependence structure of the process. 
%   X is a n-by-p array, where each row is a p-dimensional sample and n
%   is the current chain sample size.
%
%   The effective sample size of a Markov chain is the size of an i.i.d. 
%   sample with the same covariance structure as the current chain. MESS is 
%   given by MESS = n * det(LAMBDA)^(1/p) / det(SIGMA)^(1/p), where LAMBDA 
%   is the sample covariance matrix and SIGMA is an estimate of the Monte 
%   Carlo covariance matrix for the Markov chain (here obtained by batch 
%   estimation).
%
%   MESS = MULTIESS(X,SIGMA) passes optional estimate of covariance matrix
%   in Markov chain central limit theorem (CLT), as returned by a previous 
%   call to MULTIESS.
%
%   MESS = MULTIESS(X,[],B) specifies the batch size for estimation of the 
%   covariance matrix in Markov chain CLT. B can take a numeric value 
%   between 1 and n/2, or a char value between: 
%     
%   'sqroot'    B=floor(n^(1/2)) (for chains with slow mixing time; default)
%   'cuberoot'  B=floor(n^(1/3)) (for chains with fast mixing time)
%   'lESS'      pick the B that produces the lowest effective sample size
%               for a number of B ranging from n^(1/4) to n/max(20,p); this 
%               is a conservative choice
%
%   MESS = MULTIESS(X,[],B,NOFFSETS). If n is not divisible by B, SIGMA is
%   recomputed for up to NOFFSETS subsets of the data with different 
%   offsets, and the output MESS is the average over the effective sample
%   sizes obtained for different offsets (default NOFFSETS=10).
%
%   MESS = MULTIESS(X,[],B,NOFFSETS,NB) specifies the number of values of B 
%   to test when B='lESS' (default NB=200). This option is unused for other 
%   choices of B.
%
%   [MESS,SIGMA] = MULTIESS(...) returns a p-by-p covariance matrix estimate.
%
%   [MESS,SIGMA,B] = MULTIESS(...) returns the batch size used to compute 
%   the covariance SIGMA.
%
%   MULTIESS also accepts as input multiple Markov chains, either passed as 
%   a cell array (such that X{i} is the i-th Markov chain) or as a 
%   n-by-d-by-k array, where X(:,:,i) contains the samples of the i-th
%   chain. In this case, the output MESS is a 1-by-k array that reports the 
%   effective sample size of each chain. Similarly, the output SIGMA is a 
%   p-by-p-by-k array where SIGMA(:,:,i) is the covariance matrix for the 
%   i-th chain; and the output B is a 1-by-k array.
%
%   Reference: 
%   Vats, D., Flegal, J. M., & Jones, G. L. "Multivariate Output Analysis 
%   for Markov chain Monte Carlo", arXiv preprint arXiv:1512.07713 (2015).
%
%   Disclaimer: This version is still work in progress, in need of more
%   thorough testing.

% Copyright (C) 2016 Luigi Acerbi
%
% This software is distributed under the GNU General Public License 
% (version 3 or later); please refer to the file LICENSE.txt, included with 
% the software, for details.

%   Author:     Luigi Acerbi
%   Email:      luigi.acerbi@gmail.com
%   Version:    23/Jul/2016 (beta)

if nargin < 2; Sigma = []; end
if nargin < 3 || isempty(b); b = 'sqroot'; end
if nargin < 4 || isempty(Noffsets); Noffsets = 10; end
if nargin < 5; Nb = []; end

% Number of MCMC chains
if iscell(X)
    nc = numel(X);
    p = size(X{1},2);
else
    p = size(X,2);
    nc = size(X,3);
end

% SIGMA can be a cell array, but convert to standard format
if iscell(Sigma)
    temp = Sigma;
    Sigma = zeros(p,p,nc);
    for i = 1:nc; Sigma(:,:,i) = temp{i}; end
    clear temp;
end
if size(Sigma,3) == 1; Sigma = repmat(Sigma,[1 1 nc]); end

% Input check for batch size B
if ischar(b) 
    if ~any(strcmpi(b,{'sqroot','cuberoot','less'}))
        error('Unknown string for batch size. Allowed arguments are ''sqroot'', ''cuberoot'' and ''lESS''.');
    end
    if ~strcmpi(b,'less') && ~isempty(Nb)
        warning('Nonempty parameter NB will be ignored (NB is used only with ''lESS'' batch size B).');
    end
elseif isnumeric(b) && all(isfinite(b))
    if isscalar(b); b = repmat(b,[1 nc]); end
    if any(b(:) < 1) || any(b(:) > n/2)
        error('The batch size B needs to be between 1 and N/2.');
    end
    if numel(b) ~= nc
        error('The batch size B needs to be a scalar or an array of the same size as the number of chains in X.');
    end
else
    error('The batch size B needs to be either ''sqroot'', ''cuberoot'' and ''lESS'' or a number between 1 and N/2.');
end

% If no output, do a plot (only with lESS method)
plotflag = nargout == 0 && ischar(b) && strcmpi(b,'less');

% Prepare arrays
mESS = zeros(1,nc);
Sigma_out = zeros(p,p,nc);
b_out = zeros(1,nc);

% Compute multiESS separately for each chain
for i = 1:nc
    if isnumeric(b); b_in = b(i); else b_in = b; end    
    if iscell(X)
        [mESS(i),Sigma_out(:,:,i),b_out(i)] = multiESS_chain(X{i},Sigma(:,:,i),b_in,Noffsets,Nb,plotflag);
    else
        [mESS(i),Sigma_out(:,:,i),b_out(i)] = multiESS_chain(X(:,:,i),Sigma(:,:,i),b_in,Noffsets,Nb,plotflag);
    end
end

if nargout > 1; Sigma = Sigma_out; end
if nargout > 2; b = b_out; end

end

%--------------------------------------------------------------------------
function [mESS,Sigma,b] = multiESS_chain(X,Sigma,b,Noffsets,Nb,plotflag)
%MULTIESS_CHAIN Compute multiESS for a MCMC chain.

[n,p] = size(X);
if p > n
    error('More dimensions than data points, cannot compute effective sample size.');
end

if ischar(b)    
    switch lower(b)
        case 'sqroot'; b = floor(n.^(1/2));
        case 'cuberoot'; b = floor(n.^(1/3));
        case 'less'; 
            b_min = floor(n^(1/4));
            b_max = max(floor(n/max(p,20)),floor(n^(1/2)));
            if isempty(Nb); Nb = 200; end
            % Try NB log-spaced values of B from B_MIN to B_MAX
            b = unique(round(exp(linspace(log(b_min),log(b_max),Nb))));
    end
end

k = numel(b);               % Number of batch sizes
theta = mean(X,1);          % Sample mean
detLambda = det(cov(X));    % Determinant of sample covariance matrix

mESS = zeros(size(b));      % Prepare mESS
newSigma = zeros(p,p,k);    % Prepare batch Sigma matrices

% Compute mESS
for i = 1:k
    [mESS(i),newSigma(:,:,i)] = ...
        multiESS_batch(X,n,p,theta,detLambda,Sigma,b(i),Noffsets);
end

% Plot graph of mESS as a function of batch size (only for lESS)
if plotflag
    plot(b,mESS,'LineWidth',1);
    hold on;
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    box off;
    xlabel('Batch size b');
    ylabel('multiESS');
end
% Return lowest mESS
if k > 1
    [mESS,idx] = min(mESS);
    b = b(idx);
    newSigma = newSigma(:,:,idx);
end

Sigma = newSigma;

end

%--------------------------------------------------------------------------
function [mESS,Sigma] = multiESS_batch(X,n,p,theta,detLambda,Sigma,b,Noffsets)
%MULTIESS_BATCH Compute multiESS for a given batch size B.

if isempty(Sigma)
    % Compute batch estimator for SIGMA
    a = floor(n/b);
    Sigma = zeros(p,p);
    
    % Average batches over multiple offsets
    offsets = unique(round(linspace(0,n - a*b,Noffsets)));

    for j = offsets
        Y = reshape(X(j+(1:a*b),:),[b,a,p]);
        Ybar = squeeze(mean(Y,1));
        Z = bsxfun(@minus,Ybar,theta);
        for i = 1:a
            Sigma = Sigma + Z(i,:)'*Z(i,:);
        end
    end
    Sigma = Sigma*b/(a-1)/numel(offsets);    
end

mESS = n*(detLambda/det(Sigma)).^(1/p);

end
