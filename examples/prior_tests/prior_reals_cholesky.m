% prior_reals_fftma Sampling an FFT-MA type prior model
clear all,close all;clc

%% Define a CHOLESKY type Gaussian a priori model
im=1; 
prior{im}.type='cholesky';
prior{im}.type='fftma';
dx=.5;
prior{im}.x=[0:dx:10]; % X array 
prior{im}.y=[0:dx:20]; % Y array 
prior{im}.m0=10;
prior{im}.Cm='1 Sph(10,90,.25)';
prior{im}.fftma_options.constant_C=0;

x=prior{im}.x;
y=prior{im}.y;
nx=length(x);
ny=length(y);
nxy=length(x)*length(y);

randn('seed',1);



%%
figure(10);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_fftma')
try
    suptitle(sprintf('%s, Cm=%s',prior{im}.type,prior{im}.Va))
end
%% SETUP TARGET DISTRIBUTION
N=10000;
prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d_target,1,1,min(d_target),max(d_target),0);
%prior{im}.d_nscore=d_nscore;
% UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
prior{im}.o_nscore=o_nscore;
%

m=sippi_prior(prior);
sippi_plot_prior(prior,m);



return

randn('seed',1);
figure(11);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colormap(sippi_colormap(1));
suptitle(sprintf('%s using target histogram, Cm=%s',prior{im}.type,prior{im}.Va))
colorbar_shift;
print_mul('prior_reals_fftma_nscore')
return
%% MOVIE
for j=1:2
randn('seed',1);
figure(11+j);clf
set(gcf,'units','normalized','outerposition',[0 0 1 1])
prior{1}.seq_gibbs.step=0.02;

if j==1;
    fname=[mfilename,'_target','.mp4'];
    vidObj = VideoWriter(fname);
else
    prior{1}=rmfield(prior{1},'o_nscore');
    fname=[mfilename,'.mp4'];
    vidObj = VideoWriter(fname);
end
open(vidObj)
for i=1:1000;
    [m,prior]=sippi_prior(prior,m);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    colormap(sippi_colormap(1));
    caxis([8 12])
    axis image
    axis tight
    set(gca,'nextplot','replacechildren');
    %xlabel('X')
    %ylabel('Y')
     % Write each frame to the file.
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end
close(vidObj);
end

return
%%
%% trying with non-stationarity
%%

%%
V=m{im};
figure(12);clf;
[DX,DY] = gradient(V,.2,.2);
quiver(x,y,DX,DY)
axis image
return
%%
%V=zeros(size(V))+10;
d=zeros(nxy,nxy);

j=0;
for ix=1:length(x);
    for iy=1:length(y);
        j=j+1;
        t=eikonal(x,y,0,V,[x(ix),y(iy)]);
        d(j,:)=t(:)';
    end
    figure_focus(2);
    imagesc(d)
    drawnow;
end
%%
[Cm]=1-semivar_synth('.1 Cm(Nug(0) + 1 Gau(1,90,.05)',d);
for i=1:nxy
    progress_txt(i,nxy)
    for j=1:i,
        c1=Cm(i,j);
        c2=Cm(j,i);
        c=mean([c1 c2]);
        Cm(i,j)=c;
        Cm(j,i)=c;
    end
    %imagesc(Cm);drawnow;
end
        
%%
[Cm2]=precal_cov_2d(nx,ny,dx,dx,'1 Sph(10,90,.25)');
imagesc([Cm,Cm2]);axis image
%%
I=eye(nxy)*2;
figure(11);
for i=1:5
    subplot(1,5,i);
    m_sim=gaussian_simulation_cholesky(0,Cm+I,1);
    imagesc(x,y,reshape(m_sim(:,1),ny,nx))
    axis image
end
   
    

%%
return
