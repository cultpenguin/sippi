% autocorrelation: Computes the autocorrelation of a series
%
% Call
%    [d_auto,lag]=autocorrelation(d)
%
%    plot(lag,d_auto);
%
%  See also: multiESS
function [d_autoc,lag]=autocorrelation(d)

if (size(d,2)==1)
    d=d';
end
N=length(d);
dm=mean(d);

dc =conv(d-dm,fliplr(d-dm));
d_autoc=dc((N):end)/(N*var(d));
lag = 0:1:(length(d_autoc)-1);
