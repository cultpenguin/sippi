% save_wavelet: save wavelet for FDTD_fwi
%
% Call:
%    save_wavelet(data,dt,fname)
% Input:
%       - data  : 1-D data stored in the file.
%       - dt    : sample interval
%      - fname : name of the binary file (def='source.E');
% 
% By Knud Cordua (2008)                                       


function save_wavelet(data,dt,fname)

if nargin<3
  fname='source.E';
end
Fid = fopen(fname,'wb');
fwrite(Fid,length(data),'int');
fwrite(Fid,data,'double');
fwrite(Fid,dt,'double');
fclose(Fid);