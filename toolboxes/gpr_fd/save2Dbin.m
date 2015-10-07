% =========================================================================
% FUNCTION save2Dbin( fpath, fname, nx, nz, data);
% This function simply saves a 2D matrix to a binary file. It requires the
% following info: - fpath : location of the binary file
%                 - fname : name of the binary file
%                 - nx    : amount of elements in x-direction
%                 - nz    : amount of elements in z-direction (or y-)
%                 - data  : 2-D data stored in the file.
% 
% By Jacques Ernst (2006)                                        ETH Zurich
%==========================================================================
function save2Dbin(fpath,fname, nx, nz, data);

Fid = fopen(sprintf('%s%s',fpath,fname),'wb');
fwrite(Fid,nx,'int');
fwrite(Fid,nz,'int');
fwrite(Fid,data,'double');
fclose(Fid);