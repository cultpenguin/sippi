% eikonal_traveltime Computes traveltime between sources and receivers by solving the eikonal equation
%
% t=eikonal_traveltime(x,y,z,V,Sources,Receivers,iuse,type);
%
%  x,y,z : arrays defining the x, y, and z axis
%  V: velocity field, with size (length(y),length(x),length(z));
%  Sources [ndata,ndim] : Source positions
%  Receivers [ndata,ndim] : Receiver positions
%  iuse (optional): optionally only use subset of data. eg.g i_use=[1 2 4];
%  type (optional): type of eikonal solver: [1]:Fast Marching(default), [2]:FD
%
%  tmap [size(V)]: travel times computed everywhere in the velocity grid
%
%%Example (2%
%
% Example 2d traveltime compuation
%
% Example (2D):
%   x=[1:1:100];
%   y=1:1:100;
%   z=1;
%   V=ones(100,100);V(:,1:50)=2;
%   S=[50 50 1;50 50 1];
%   R=[90 90 1; 90 80 1];
%   t=eikonal_traveltime(x,y,z,V,S,R)
%
% Example (3D):
%   nx=50;ny=50;nz=50;
%   x=1:1:nx;
%   y=1:1:ny;
%   z=1:1:nz;
%   V=ones(ny,nx,nz);V(:,1:50,:)=2;
%   S=[10 10 1;10 10 1;10 9 1];
%   R=[40 40 40; 40 39 40; 40 40 40];
%   t=eikonal_traveltime(x,y,z,V,S,R)
%
%
% See also eikonal

function t=eikonal_traveltime(x,y,z,V,S,R,iuse,type);

if nargin<7,
    ns=size(S,1);
    iuse=1:1:ns;
end
if nargin<8, type=1;end
ns=length(iuse);


if size(S,2)<3;
    if length(z)==1;
        S(:,3)=z;
    end
end
if size(R,2)<3;
    if length(z)==1;
        R(:,3)=z;
    end
end


S=S(iuse,:);
R=R(iuse,:);


Su=unique(S,'rows');

tmap=eikonal(x,y,z,V,Su,type);
tmap=squeeze(tmap); %% Output 'type=1' is 4D, and for 'type=2' 3D

[xx,yy,zz]=meshgrid(x,y,z);


t=zeros(ns,1);
for i=1:size(Su,1);
    
    if length(z)==1
        % 2D
        ir=find( (S(:,1)==Su(i,1)) & (S(:,2)==Su(i,2))  & (S(:,3)==Su(i,3)) );
        t(ir)=interp2(xx,yy,tmap(:,:,i),R(ir,1),R(ir,2));
    else
        % 3D
        ir=find( (S(:,1)==Su(i,1)) & (S(:,2)==Su(i,2))  & (S(:,3)==Su(i,3)) );
        t(ir)=interp3(xx,yy,zz,tmap(:,:,:,i),R(ir,1),R(ir,2),R(ir,3));
        ii=find(isnan(t(ir)));
        if ~isempty(ii)
            keyboard
        end
        
    end
end

