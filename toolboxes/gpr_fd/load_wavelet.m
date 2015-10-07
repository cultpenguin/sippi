% load_wavelet: load wavelet for FDTD_fwi
%
% Call: [data,dt]=load_wavelet(fname);
%
% Input:
%    - fname [def='source.E']: name (and optionally full path) of the binary file
% Output
%    - data  : 1-D data stored in the file.
%    - dt    : sample interval

function [data,dt]=load_wavelet(fname)

if nargin==0
  fname='source.E';
end

Fid = fopen(fname,'rb');
n = fread(Fid,1,'int');
data = fread(Fid,'double');
dt=data(end);
data=data(1:(end-1));
fclose(Fid);