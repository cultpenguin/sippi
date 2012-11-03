% eikonal Traveltime computation by solving the eikonal equation
% 
% tmap=eikonal(x,y,z,V,Sources,type);
%
%  x,y,z : arrays defining the x, y, and z axis
%  V: velocity field, with size (length(y),length(x),length(z));
%  Sources [ndata,ndim] : Source positions
%  type (optional): type of eikonal solver: [1]:Fast Marching(default), [2]:FD
%
%  tmap [size(V)]: travel times computed everywhere in the velocity grid
%
%%Example (2D):
%   x=[1:1:100];
%   y=1:1:100;
%   z=1;
%   V=ones(100,100);V(:,1:50)=2;
%   Sources = [10 50;75 50];
%   t=eikonal(x,y,z,V,Sources);
%   subplot(1,2,1);imagesc(x,y,t(:,:,1,1));axis image;colorbar
%   subplot(1,2,2);imagesc(x,y,t(:,:,1,2));axis image;colorbar
%
% See also eikonal_traveltime
%
function tmap=eikonal(x,y,z,V,Sources,type);

if nargin<6, type=1; end
if nargin<5, 
    Sources=[x(1),y(1),z(1)];
end

ny=size(V,1);
nx=size(V,2);
nz=size(V,3);
ns=size(Sources,1);
dx=x(2)-x(1);

if size(Sources,2)<3;
    if nz==1;
        Sources(:,3)=z;
    end
end
    

if type==1
    % FAST MARCHING
    
    SourcesUnique=unique(Sources,'rows');
    nu=size(SourcesUnique,1);
    if nu<ns
        % ONLY COMPUTE 
        tmapU=eikonal(x,y,z,V,SourcesUnique,type);
        tmap=zeros(ny,nx,nz,ns);
        for j=1:ns
            iu=find( (SourcesUnique(:,1)==Sources(j,1)) & (SourcesUnique(:,2)==Sources(j,2)) & (SourcesUnique(:,3)==Sources(j,3)));
            tmap(:,:,j)=tmapU(:,:,iu);
        end
        return
    end
    
    % CONVERT TO INDEX LOCATIONS
    Sources_ind=Sources.*0;
    for j=1:size(Sources,2);
       if j==1; Sources_ind(:,j)=((Sources(:,j)-x(1))./dx)+1;end
       if j==2; Sources_ind(:,j)=((Sources(:,j)-y(1))./dx)+1;end
       if j==3; Sources_ind(:,j)=((Sources(:,j)-z(1))./dx)+1;end
    end
    
    %
    tmap=zeros(ny,nx,nz,ns);
    
    for is=1:size(Sources_ind,1);
        %disp(is)
        % SWAP X and Y
        S=Sources_ind(is,:);
        S(1)=Sources_ind(is,2);
        S(2)=Sources_ind(is,1);
        %
        [tmap(:,:,:,is)]=msfm(V, round(S'), true, true).*dx;
    end
   
    
elseif type==2
    % NFD / EIKONAL
    if nz==1;
        % 2D     
        tmap = fast_fd_2d(x,y,V,Sources(:,1:2));
    else
        % 3D
        disp('3D not implemented for NFD')
        tmap=[];
    end
end    
