% prob_boundary_1d(D,boundary_width)
%
% IN:
%   D [nx,n_reals] : matrix with n_reals realizations of length nx
%   boundary_width [real]: minimum boundary change
% OUT: 
%   prob_boundary [nx-1,1]: probability of change>boundary_width
%
% Call :
%   prob_boundary=prob_boundary_1d(D,boundary_width)
%
%
function prob_boundary=prob_boundary_1d(D,boundary_width)

if nargin<1
    help(mfilename);
    n_boundary=[];
end
[nx,nt]=size(D);


if nargin<2
    boundary_width=[max(D(:))-min(D(:))]/10;
end

prob_boundary=sum([abs(diff(D))>boundary_width]')'/nt;

