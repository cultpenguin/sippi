% eikonal_traveltime : computes traveltime by solving the eikonal equation
%
% t=eikonal_traveltime(x,y,z,V,Sources,Receivers,iuse,type);
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

