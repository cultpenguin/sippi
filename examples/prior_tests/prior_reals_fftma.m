% prior_reals_fftma Sampling an FFT-MA type prior model
clear all,

%% Define an FFTMA type a priori model
im=1; 
prior{im}.type='FFTMA';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
prior{im}.Va='1 Sph(10,90,.25)';
prior{im}.fftma_options.constant_C=0;

prior=sippi_prior_init(prior);

randn('seed',1);
%%
figure(10);clf
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
colorbar_shift;
print_mul('prior_reals_fftma')
suptitle(sprintf('FFT-MA, Cm=%s',prior{im}.Va))

%% SETUP TARGET DISTRIBUTION
N=10000;
prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
[d_nscore,o_nscore]=nscore(d_target,1,1,min(d_target),max(d_target),0);
% UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
prior{im}.o_nscore=o_nscore;
prior=sippi_prior_init(prior);
%
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
suptitle(sprintf('FFT-MA using target histogram, Cm=%s',prior{im}.Va))
colorbar_shift;
print_mul('prior_reals_fftma_nscore')

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
