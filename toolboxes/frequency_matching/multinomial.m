function logP = multinomial(Htest, Htrain,prior)

% log probability using the dirichlet distribution
% It becomes the multinomial distribution when prior =0 
% H has to be coloumn vectors.
%
% KC 2016
%
if nargin<3
    prior=0;
end

index = Htrain(:) > 0 | Htest(:) > 0;
Htrain = Htrain(index);
Htest = Htest(index);

Hprior = prior*ones(size(Htrain));

Nprior = sum(Hprior);
Ntrain = sum(Htrain);
Ntest = sum(Htest);
Nopt = Nprior + Ntrain;

Hopt = (Hprior + Htrain);

tmp1 = gammaln(Htest + 1);
tmp3 = log( Hopt / Nopt );

logP = gammaln(Ntest + 1) - sum(tmp1) + Htest' * tmp3;

%P=exp(P);

%Pprob = factorial(Ntest) / prod(factorial(Htest(:))) * prod((Hopt(:)/Nopt).^Htest(:))

% Ptest = gammaln(Ntest + 1) - sum(gammaln(Htest(:)+1)) + sum(Htest(:).*log(Htrain(:)/Ntest))

