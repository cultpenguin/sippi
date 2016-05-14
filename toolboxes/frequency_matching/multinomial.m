function [loglik,lik] = multinomial(H,Hti,Hprior,type)
if nargin==0;
    H=[2 3 5 2 1]*12;
    Hti=[2 3 5 2 1]*22;
    Hprior=[1 1 1 1 1].*0;
end
if nargin<3;
    Hprior=1;
end

if length(Hprior)==1
    Hprior=ones(size(H)).*Hprior;
end

if nargin<4
    type=1;
end

N=sum(H);
Nti=sum(Hti);
Nprior=sum(Hprior);
Nb=length(H);




% Cordua et al 2014



if type==1;
    %% LOG-probability
    K_top = sum(log(1:N));
    for i=1:Nb
        K_bas(i) = sum(log(1:H(i)));
    end
    K= K_top - sum(K_bas);
    
    
    for i=1:Nb
        Li(i)= H(i)*log( (Hti(i) + Hprior(i))/(Nti+Nprior) );
    end
    L=NaN;
    L=sum(Li);
    loglik=K+L;
    lik=exp(loglik);
elseif type==2
    %% probability
    P_top = factorial(N);
    for i=1:Nb
    P_bas(i) = factorial(H(i));   
    end
    P = P_top./prod(P_bas);
    
    for i=1:Nb
        Qi(i)= ( (Hti(i) + Hprior(i))/(Nti+Nprior) )^H(i);
    end
    Q=prod(Qi);
    
    lik=P*Q;
    loglik=log(lik);
end

%%
%[log(P),K]
%[log(Q),L]
%%keyboard
%[log(lik),loglik]