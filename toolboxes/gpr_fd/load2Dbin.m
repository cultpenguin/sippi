% =========================================================================
% FUNCTION [ nx, nz, data] = load2Dbin( fpath, fname);
% This function simply loads a 2D matrix from a binary file. It requires
% the following info: - fpath : location of the binary file
%                     - fname : name of the binary file
% Remark: if only one input argument is given, the function assumes that it
% contains info about path AND filename!!
%
% And it returns:     - nx    : amount of elements in x-direction
%                     - nz    : amount of elements in z-direction (or y-)
%                     - data  : 2-D data stored in the file.
% 
% By Jacques Ernst (2006)                                        ETH Zurich
%==========================================================================
function [ nx, nz, data] = load2Dbin( fpath, fname);

if (nargin==1)
    tmpp = sprintf('%s',fpath); txt = 'Noname';
else
    tmpp = sprintf('%s\\%s',fpath,fname); txt = fname;
end;

Fid = fopen(sprintf('%s',tmpp),'rb');
if (Fid==-1)
    fclose all;
    data = []; nx = []; nz = [];
    disp(sprintf('WARNING: %s-File not found...',txt));
else
    nx = fread(Fid,1,'int');
    nz = fread(Fid,1,'int');
    data = fread(Fid,[nz nx],'double');
    fclose(Fid);
end;