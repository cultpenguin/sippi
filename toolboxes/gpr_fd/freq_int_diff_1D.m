function odata=freq_int_diff_1D(idata,order,dt)

% Call: [odata opwr]=freq_int_diff_1D(idata,order,dt); 
% Outputs:
% * "odata" is the output data
% * "opwr" is the log power spectrum of the output data
% Inputs:
% * "idata" is the input data given as a matrix
% * "order" A positive value N results in a differentialtion of order N. A
%    negative value N results in an integration of order N. (e.g. for N=1/2 the array is half differentiated. N = -2 results in a double integration)
% * "dt" is the sample interval
% 
% The code will extend the input signal with zeroes to the next power of two
% 
% For finite-difference calculations; if a e.g. Ricker pulse has to propagte in
% the grid then a half-integrated (order = -0.5) Ricker pulse has to be
% used as source pulse.
%
% Knud S. Cordua 2010

if order<0
    if abs(order)==round(abs(order))
        % IN CASE order=-1, -2, -3, -4, ..
        n_diff=0;
        n_int=abs(order);
    else
        n_diff=1-mod(abs(order),1);
        n_int=ceil(abs(order));
    end
else 
    n_diff=order;
end

N=2^nextpow2(length(idata));
%N=length(idata);

% Differentiation
wf=fft(idata,N);
wf=fftshift(wf);
f=disc_freq(N,dt);
ii=sqrt(-1);
uu=2*pi*f;
wf=((ii.*uu).^n_diff).*wf;
wf=fftshift(wf);
odata=real(ifft(wf));

% Integration
if order<0
    for i=1:n_int;
    odata=cumsum(odata)*dt;
    end
end


