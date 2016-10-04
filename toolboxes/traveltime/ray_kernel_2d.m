function A=ray_kernel_2d(antpos,Nz,Nx,cell,r,norm)

% Calculates a ray based sensitivity kernel:
%
% Call: A=calc_G(ant_pos,Nz,Nx,dim,r,norm);
% The input r and norm are optional
% 
% "nz" and "nx" are the number of model parameters in depth and distance
% respectively. "cell" is the size of the cells given in meters.
% "ant_pos" is a d x 4 matrix, where d is the number of data. The 4 columns 
% represents the x1 z1 x2 z2 coordinates of the transmitter and receiver 
% position in meters
% "r" is a refinement factor. If r=2 the kernel will be calculated with
% cells with sides half as large as the input. This is usedfull when odd
% positions with respect to the grid is needed.
% "norm" indicates if the sensitivity kernel is normalised (=1) or not (=0)
%
% Originally developed by Bo Holm Jacobsen

if nargin<6, norm=0; end
if nargin<5, r=1;norm=0; end

dim=cell/r;
Nz=Nz*r;
Nx=Nx*r;
shotrec=antpos./dim;
if r~=1, disp([mfilename,': Remember to move the antennae positions to the midle of the refined cells']),
    disp(sprintf('%s: %s %d',mfilename,'The refined cell size is',dim))
end

fig_num=gcf;

xsps=shotrec(:,3); zsps=shotrec(:,4);
xrs=shotrec(:,1);  zrs=shotrec(:,2);

x0 = 0; dx=1;  xNx=x0+Nx*dx;
z0 = 0; dz=1;  zNz=z0+Nz*dz;

show_rays=0;


if show_rays==1,
figure
plot(0,0,'k'), xlabel('Distance (Cell no.)','fontweight','demi'), ylabel('Depth (Cell no.)','fontweight','demi');hold on
end

%  for nx=0:Nx
%    plot([x0+nx*dx,x0+nx*dx],-[z0,zNz],'k');
%  end
%  for nz=0:Nz
%    plot([x0,xNx],-[z0+nz*dz,z0+nz*dz],'k');
%  end
 % hstar=line(xsps,-zsps);
 %plot(xsps,-zsps,'r*');
 %plot(xrs,-zrs,'go');
Nray=length(xsps);

%******defining the geometry of the rays***********
A = ones(Nray,Nx*Nz);
Arow = zeros(Nz,Nx);

xs=x0; xr=x0; zs=z0; zr=z0; % just to allow plotting of rays

%************start of loop************
for iray=1:Nray

if show_rays==1,
  hpre=line([xs,xr],-[zs,zr]); % over-printing previous ray in blue
  %set(hpre,'linewidth',0.0000001,'color','r','linestyle','--')
  set(hpre,'linewidth',1.5,'color','r','linestyle','--');
end

xs=xsps(iray);
zs=zsps(iray);
xr=xrs(iray);
zr=zrs(iray);

if show_rays==1,
  hcur=line([xs,xr],-[zs,zr]); % showing current ray in green
  %set(hcur,'linewidth',0.0000001,'color','r','linestyle','--');
  set(hcur,'linewidth',1.5,'color','r','linestyle','--');
end

if xs==xr, % vertical ray; we tilt it slightly
% reprogram to allow for very different magnitude of zs and zr !!
  xs=(xs+100000*pi*eps)*(1+ 1234*eps);
  xr=(xr-100000*pi*eps)*(1- 1234*eps);   
end

% if zs==zr, % horizontal ray; we tilt it slightly
%   zs=(zs+100000*pi*eps)*(1+ 1234*eps);
%   zr=(zr-100000*pi*eps)*(1- 1234*eps);   
% 
% end

zs_tmp=zs;
zr_tmp=zr;

if zs==zr
    zs=zs+10^-10;
    zr=zr-10^-10;
end
    

if zr<zs, x1=xr; z1=zr; x2=xs; z2=zs;    
  else    x1=xs; z1=zs; x2=xr; z2=zr; 
end

% solve for a and b for line through (x1,z1) and (x2,z2).
ab = [x1,1;x2,1]\[z1;z2]; 
a=ab(1); b=ab(2);

%.... projection external source and/or receiver along ray onto grid boundary
x_x0  = x0;        z_x0 =a*x0 +b; % coordinates of ray intersection with x=x0
x_z0  = (z0-b)/a;  z_z0 =z0 ;     % coordinates of ray intersection with z=z0
x_xNx = xNx;       z_xNx=a*xNx+b; % coordinates of ray intersection with x=x0
x_zNz = (zNz-b)/a; z_zNz=zNz;     % coordinates of ray intersection with z=z0

if a>0 
  x1eff = max([x1;x_x0 ;x_z0 ]);
  x2eff = min([x2;x_xNx;x_zNz]); 
else
  x1eff = min([x1;x_xNx;x_z0 ]);
  x2eff = max([x2;x_x0 ;x_zNz]);
end

z1eff = a*x1eff + b;
z2eff = a*x2eff + b;

% ................. effective points are moved infinitesimally inside grid
z1eff = z1eff + eps*(1+abs(z1eff)); % good except for atomic distances
z2eff = z2eff - eps*(1+abs(z2eff));
% zmove = eps * max(abs(x0),abs(xNx) * abs(a);
z1eff = z1eff + eps*Nx*dx*abs(a);
z2eff = z2eff - 10*eps*dx*abs(a);
x1eff = (z1eff-b)/a;
x2eff = (z2eff-b)/a;

ix1 = ceil((x1eff-x0)/dx);  
%ix1  = max(1,min(Nx,ix1)); % index of first column entered by ray; 0 and Nx suppressed.
    
ix2 = ceil((x2eff-x0)/dx);   
%ix2  = max(1,min(Nx,ix2)); % index of last  coulmn entered by ray; 0 and Nx suppressed.  

zout = z1eff; 
izout = ceil((zout-z0)/dz);
xout = x1eff;
Arow = zeros(Nz,Nx); 

for ix=ix1:sign(a):ix2
  % deal with ray segments between x=x0+(ix-1)*dx and x0+ix*dx
  xin = xout;
  zin = zout;
  if a>0
    xout = min(x2eff,x0+ix*dx);
  else
    xout = max(x2eff,x0+(ix-1)*dx);
  end
    
  zout = a*xout + b;
  izin  = izout;
  izout = ceil((zout-z0)/dz);
  
  if izout>=izin, % at least one cell was hit (relevant for certain rays outside grid)
    zz1 = min(zout,z0+izin*dz) - zin; % z-increment in first cell
    xx1 = zz1/a;                      % x-increment in first cell
    L1  = sqrt(zz1*zz1 + xx1*xx1);    % ray length in first cell
    Arow(izin,ix) = L1;
  end

  if izout>izin, % more than one cell in column ix is hit
    zzm = zout-(z0+(izout-1)*dz);
    xxm = zzm/a;
    Lm  = sqrt(zzm*zzm + xxm*xxm);
    Arow(izout,ix) = Lm; 
  end % if izout>izin

  if izout>(izin+1), % hitting cells between first and last hit cell in column.
    zz=dz;
    xx=zz/a;
    L = sqrt(zz*zz+xx*xx);
    if izout==izin+2
      Arow(izin+1,ix) = L;
    else
      Arow(izin+1:izout-1,ix) = L*ones(size(izin+1:izout-1))'; 
    end
  end % if izout>(izin+1)
end % for ix=ix1:ix2


if zs_tmp==zr_tmp
    AA=Arow(:);
    for i=1:length(AA)
        if AA(i)>0
            AA(i)=1;
        end
    end
    A(iray,:) =AA';
else
A(iray,:) = Arow(:)';
end

end % for iray=1:Nray

A=A.*dim;
if show_rays==1,
axis image

for nx=0:Nx
   plot([x0+nx*dx,x0+nx*dx],-[z0,zNz],'k');
end
for nz=0:Nz
   plot([x0,xNx],-[z0+nz*dz,z0+nz*dz],'k');
end
end
if norm==1
    for i=1:length(A(:,1))
        A(i,:)=A(i,:)./sum(A(i,:));
    end
end

if r~=1
   for k=1:length(A(:,1))
        A1=reshape(A(k,:),Nz,Nx);
        for i=1:Nz/r
            for j=1:Nx/r
                A2=A1((((i-1)*r+1):i*r),((j-1)*r+1):j*r);
                A3(i,j)=sum(A2(:));
            end
        end
        A_red(k,:)=A3(:);
   end
   A=A_red;
end


%close figure 1