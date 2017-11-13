% autocorrelation
%
%  See also: multiESS
function d_autoc=autocorrelation(d)

if (size(d,2)==1)
    d=d';
end
N=length(d);
dm=mean(d);

dc =conv(d-dm,fliplr(d-dm));
d_autoc=dc((N+1):end)/(N*var(d));

