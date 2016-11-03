
function [data time_vector]=image_snapshots(sample_rate,Nplot,Ntrn,Nrow,gain,dx,dt,field)

%-------------------- Make plot of snapshots ------------------------------
%
% Call: data=image_snapshots(sample_rate,Nplot,Ntrn,Nrow,gain,dx,dt,field);
%
% * sample_rate: The rate of time-steps at which snapshots are written to
%   the disk during the FDTD forward modelling
% * Nplot: The number of snapshots to be plotted starting from 1, or the interval of
%   snapshots to be plotted given as a vector.
% * Ntrn: The transmitter number (starting with 0) from which the snapshots
%   are plotted
% * Nrow: Numer of plots/row
% * gain: A gain factor used to amplify the snapshots as function of the
%   snapshot number. The amplification is obtained by scaling the initial
%   color axis with the factor g=(i*(gain+1))^-2. For a gain factor of 0
%   now amplifcation is applied to the snapshot.
% * dx: The cell size (only needed to obtain the correct axis in the images)
% * dt: Used to determine the time of the snapshots
% * field: The direction of the eletrical field 1: Ez, 2: Ex 
%
% Knud Cordua (2008)
%
%--------------------------------------------------------------------------

do_plot=1;

if nargin==3
    Nrow=8;
    gain=1;
end

if nargin==4
    gain=1;
end

init_caxis=[-1 1];

if length(Nplot)==1
    Np(1)=1;
    Np(2)=Nplot;
    Nplot=Np;
end


N_vert=ceil((max(Nplot)-min(Nplot)+1)/Nrow);

c=0;
%colormap(bone)
time_vector=[];
%figure
for i=min(Nplot):max(Nplot)
    c=c+1;

    time=i*sample_rate;
    fname=timetag(time,Ntrn,field);
    try
        [data(:,:,c),nx,nz]=load2Dbin(fname);
    catch
        disp(sprintf('%s %i %s','Warning: only',i-1','plots are available (N_plots too large)'))
        return
    end

    if do_plot
    subplot(N_vert,Nrow,i-min(Nplot)+1);
    
    g=(i*(gain+1))^-2;
    if nargin<6
        imagesc(data(:,:,c)),axis image,caxis(init_caxis*g)
        xlabel('Distance (cells)'),ylabel('Depth (cells)')
    else
        imagesc(0:dx*nx,0:dx*nz,data(:,:,c)),axis image,caxis(init_caxis*g)
        xlabel('Distance (m)','fontsize',16),ylabel('Depth (m)','fontsize',16)
    end
    if nargin<7
        title(sprintf('%1.0f %s %1.1f',time,'time-steps, Gain=',1/g))
    else
        title(sprintf('%1.0f %s %1.1f',dt*time*10^9,'ns'))
    end
    end
    time_vector=[time_vector dt*time*10^9];
end


function out=timetag(tag,Ntrn,field)
if field==1
    out=sprintf('forward_EzS%05d_T%02d.dat',tag,Ntrn);
else
    out=sprintf('forward_ExS%05d_T%02d.dat',tag,Ntrn);
end

function [data,nx,nz]=load2Dbin(fname)
fid = fopen(fname,'rb');
nx = fread(fid,1,'int');
nz = fread(fid,1,'int');
data = fread(fid,[nz nx],'double');
fclose(fid);
