function p = normpdf(x,mu,cov,varargin)

N = length(x);

if strcmp(varargin,'log')
    p = -(N/2)*log(2*pi) - (1/2)*logdet(cov) ...
        - (1/2) * (x-mu)' * inv(cov) * (x-mu);
else 
    p = ( 1/sqrt( ((2*pi)^N) * det(cov) ) ) *...
        exp( -.5 * (x-mu)' * inv(cov) * (x-mu));
end

end
