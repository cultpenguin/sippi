% tomography_kernel Computes the sensitivity kernel for a wave traveling from S to R.
%
% CALL :
%    [K,RAY,Gk,Gray,timeS,timeR,raypath]=tomography_kernel(Vel,x,y,z,S,R,T,alpha,Knorm);
%
% IN :
%    Vel [ny,nx] : Velocity field
%    x [1:nx] :
%    y [1:ny] :
%    z [1:nz] :
%    S [1,3] : Location of Source
%    R [1,3] : Location of Receiver
%    T : Donminant period
%    alpha: controls exponential decay away ray path
%    Knorm [1] : normaliztion of K [0]:none, K:[1]:vertical
%
% OUT :
%    K : Sensitivity kernel
%    R : Ray sensitivity kernel (High Frequency approx)
%    timeS : travel computed form Source
%    timeR : travel computed form Receiver
%    raypath [nraydata,ndim] : the center of the raypath
%
% The sensitivity is the length travelled in each cell.
%
%
%
function [K,RAY,Gk,Gray,tS,tR,raypath_mat,raylength_mat]=tomography_kernel(Vel,x,y,z,S,R,T,alpha,Knorm,doPlot);
if nargin<7, T=1; end
if nargin<8, alpha=1; end
if nargin<9,
    Knorm=1;
end
if nargin<10,
    doPlot=0;
end

nx=length(x);
ny=length(y);
nz=length(z);

if nz==1;is3d=0;else;is3d=1;end
    

if length(Vel(:))==1;
    Vel=ones(length(y),length(x),length(z)).*Vel;
end
x0=x(1);
y0=y(1);
z0=z(1);
dx=x(2)-x(1);

ns=max([size(S,1) size(R,1)]);
dx=x(2)-x(1);
dy=y(2)-y(1);
d1=(dx+dy)/2;

itype=1; % MSFM
%itype=2; % NFD / FAST
tS=eikonal(x,y,z,Vel,S,itype);
tR=eikonal(x,y,z,Vel,R,itype);

if (size(tS,3)==1)*(size(tR,3)>1)
    disp('!!!')
    ttS=tR.*0;
    for i=1:size(tR,3)
        ttS(:,:,i)=tS;
    end
    tS=ttS;
    S=repmat(S,[ns 1]);
end
if exist('fast_fd_clean','file')
    fast_fd_clean;
end

SR_dist_linear=sqrt(sum((S-R).^2')');


dt=tS+tR;
K=zeros(size(dt));
RAY=zeros(size(dt));
str_options = [.01 50000]; % CONTROL STREAM2 BELOW!!
[xx,yy]=meshgrid(x,y);
for is=1:ns
    
    if ((is/10)==round(is/10))
        progress_txt(is,ns,mfilename);
    end

    % geometrical spreading type
    % spread_type=0; % PLANE
    spread_type=1; % CYLINDRICAL
    % spread_type=2; % SPHERICAL
    
    if is3d
        mt=min(min(min(dt(:,:,:,is))));
        dt(:,:,:,is)=dt(:,:,:,is)-mt;
        % GEOMETRICAL SPREADING
        aS=tS(:,:,:,is);aS(find(aS==0))=d1;
        aR=tR(:,:,:,is);aR(find(aR==0))=d1;
        aR=spherical_spreading(aR,spread_type);
        aS=spherical_spreading(aS,spread_type);
        % CALCULATE KERNEL
        %K(:,:,is)=munk_fresnel_3d(T,dt(:,:,is),alpha,aS,aR);
        K(:,:,:,is)=munk_fresnel_3d(T,dt(:,:,:,is),alpha);
    
    
    else
        mt=min(min(dt(:,:,is)));
        dt(:,:,is)=dt(:,:,is)-mt;
        % GEOMETRICAL SPREADING
        aS=tS(:,:,is);aS(find(aS==0))=d1;
        aR=tR(:,:,is);aR(find(aR==0))=d1;
        aR=spherical_spreading(aR,spread_type);
        aS=spherical_spreading(aS,spread_type);
        % CALCULATE KERNEL
        %K(:,:,is)=munk_fresnel_2d(T,dt(:,:,is),alpha,aS,aR);
        K(:,:,is)=munk_fresnel_2d(T,dt(:,:,is),alpha);
    
    end
    
    % NOW FIND FIRST ARRIVAL AND RAYLENGTH
    [U,V]=gradient(tS(:,:,is));
    start_point=R(is,:);
    raypath = stream2(xx,yy,-U,-V,start_point(1),start_point(2),str_options);
    raypath=raypath{1};
    
    % MAKE CHECK THAT THAT THE RAY IS ACTUALLY TRACED WITHIN A SMALL
    % DISTANCE OF THE SOURCE
    end_point=S(is,:);
    
    dist_to_source=sqrt(sum( (raypath(end,:)-S(is,:)).^2));
    if (dist_to_source>dx)||(dist_to_source>dy);
        if (dist_to_source<2*dx)&&(dist_to_source<2*dy);
            sippi_verbose(sprintf('%s: Unstable ratracing: S=[%g,%g], R=[%g,%g] (%g,%g)',mfilename,S(is,1),S(is,2),start_point(1),start_point(2)),-1)
            sippi_verbose(sprintf('%s: CHECK OPTIONS FOR ''stream2''',mfilename),-1)
        else
            sippi_verbose(sprintf('%s: Unable to trace the ray to the source point for S=[%g,%g], R=[%g,%g] (%g,%g)',mfilename,S(is,1),S(is,2),start_point(1),start_point(2)),-1)
            sippi_verbose(sprintf('%s: CHECK OPTIONS FOR ''stream2''',mfilename),-1)
            sippi_verbose(sprintf('%s: RESULTS MAY NOT BE USEFUL!!!''',mfilename),-1)
            
        end
    end
    
    if isempty(raypath)
        sippi_verbose(sprintf('%s: Unable to compute raypath from start point (%g,%g)',mfilename,start_point(1),start_point(2)),-1)
    end
    
    % GET RID OF DATA CLOSE TO SOURCE (DIST<DX)
    r2=raypath;r2(:,1)=r2(:,1)-S(is,1);r2(:,2)=r2(:,2)-S(is,2);
    distS=sqrt(r2(:,1).^2+r2(:,2).^2);
    ClosePoints=find(distS<dx);
    
    % NEXT LINES OBSOLETE
    ic=0;
    while (isempty(ClosePoints))
        ic=ic+1;
        ClosePoints=find(distS<(ic*dx/10));
        if ic==20;break;end
    end
    
    if isempty(ClosePoints)
        igood=1:1:length(distS);
        % something went wrong
    else
            igood=1:1:ClosePoints(1);
    end
    
    raypath=[raypath(igood,:)];
    % ADD Source Locations as last point;
    raypath=[raypath;S(is,1:2)];
    raylength=sum(sqrt(diff(raypath(:,1)).^2+diff(raypath(:,2)).^2));
    
    force_linear_raylength=0; % OK IN CASE OF CONST VEL FIELD
    if force_linear_raylength==1;
        % force linear raylength
        raylength=SR_dist_linear(is);
    end
    
    
    raypath_mat{is}=raypath;
    raylength_mat(is)=raylength;
    
    ix=ceil((raypath(:,1)-(x0-dx/2))./dx);
    iy=ceil((raypath(:,2)-(y0-dx/2))./dx);
    
    ix(find(ix<1))=1;
    iy(find(iy<1))=1;
    
    
    for j=1:length(ix)
        RAY(iy(j),ix(j),is)=RAY(iy(j),ix(j),is)+1;
    end
    
    sk=K(:,:,is);sk=sum(sk(:));
    K(:,:,is)=raylength_mat(is).*K(:,:,is)./sk;
    
end

% NORMALIZE RAY
for is=1:ns
    r=RAY(:,:,is);
    RAY(:,:,is)=r.*raylength_mat(is)./sum(r(:));
end


if Knorm==1
    % Vertical normalization of Fresenel Kernel
    % SIMPLE normalization in case of cross borehole
    % inversion between two vertical borehole
    % i.e. when wave travel horizontally.
    for is=1:ns
        single_ray=RAY(:,:,is);
        sRAY=sum(single_ray);
        for i=1:size(single_ray,2);
            sk=sum(K(:,i,is));
            if sk>0
                K(:,i,is)=sRAY(i).*K(:,i,is)./sk;
            end
        end
    end
end
% REPORT
Gray=zeros(ns,length(x)*length(y));
Gk=zeros(ns,length(x)*length(y));
for is=1:ns
    g=RAY(:,:,is);
    %g=g./sum(g(:));
    Gray(is,:)=g(:);
    
    gk=K(:,:,is);
    %    gk=gk./sum(gk(:));
    Gk(is,:)=gk(:);
end


if doPlot>1;
    figure;
    ip=size(K,3);
    subplot(2,3,1)
    imagesc(x,y,Vel);axis image;title('Velocity model')
    subplot(2,3,2)
    imagesc(x,y,tS(:,:,ip));axis image;title('t_{source}')
    subplot(2,3,3)
    imagesc(x,y,tR(:,:,ip));axis image;title('t_{receiver}')
    
    subplot(2,3,4)
    imagesc(x,y,K(:,:,ip));axis image;title('Fresnel kernel')
    hold on
    plot(S(ip,1),S(ip,2),'r*')
    plot(R(ip,1),R(ip,2),'ro')
    plot(raypath(:,1),raypath(:,2),'w*','Markersize',2)
    colorbar
    hold off
    
    subplot(2,3,5)
    imagesc(x,y,RAY(:,:,ip));axis image
    hold on
    plot(S(ip,1),S(ip,2),'r*')
    plot(R(ip,1),R(ip,2),'ro')
    hold off
    colorbar
    title('Ray kernel')
    drawnow;
end

if doPlot>2;
    for i=1:ns;imagesc(dt(:,:,i));axis image;drawnow;end
    for i=1:ns;imagesc(K(:,:,i));axis image;drawnow;end
    for i=1:ns;imagesc(RAY(:,:,i));axis image;drawnow;end
end


