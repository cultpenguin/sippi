% prior_reals_visim Sampling a VISIM type prior model

clear all,
im=1; 
prior{im}.type='VISIM';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
%% UPDATE : MAKE SURE THIS WORKS FOR SPECIFIC COV.
prior{im}.Va='1 Sph(10,90,.25)'; 

prior=sippi_prior_init(prior);

[m,prior]=sippi_prior(prior);
prior{1}.V.Va.a_hmax=10;
prior{1}.V.Va.a_hmin=2.5;
prior{1}.V.Va.ang1=90;
prior{1}.V.Va.it=1;

randn('seed',1);rand('seed',1);
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
try;load cmap;colormap(cmap);end
colorbar_shift
print_mul('prior_reals_visim')



%%
N=10000;
prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
%prior{im}.d_target=d_target;
[d_nscore,o_nscore]=nscore(d_target,1,1,min(d_target),max(d_target),0);
prior{im}.o_nscore=o_nscore;
prior{im}.d_target=d_target;
prior=sippi_prior_init(prior);
%
randn('seed',1);rand('seed',1);
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
try;load cmap;colormap(cmap);end
colorbar_shift
print_mul('prior_reals_visim_target')

%% SEQ GIBBS TYPE 1
try;prior{im}=rmfield(prior{im},'d_target');end
randn('seed',2);rand('seed',2);
prior{im}.seq_gibbs.type=1;
prior{im}.seq_gibbs.step=[4];

%m_old=m;
%[m,prior]=sippi_prior(prior,m);
%imagesc(m{1}-m_old{1});colorbar;axis image;caxis([-1 1].*1e-5)
[m,prior]=sippi_prior(prior);

figure(12);clf
for i=1:5;
    m_old=m;
    [m,prior]=sippi_prior(prior,m);
    subplot(2,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
     ax=axis;
    subplot(2,5,5+i);
    d_diff=m{1}.*0;
    ii=find(m_old{1}~=m{1});
    plot(prior{1}.xx(ii),prior{1}.yy(ii),'k.');
    set(gca,'ydir','revers')
    axis(ax)
    
end
try;load cmap;colormap(cmap);end
subplot(2,5,5);colorbar_shift;
print_mul('prior_reals_visim_seqgibbs_type1');


%% SEQ GIBBS TYPE 2
try;prior{im}=rmfield(prior{im},'d_target');end
randn('seed',1);rand('seed',1);
prior{im}.seq_gibbs.type=2;
prior{im}.seq_gibbs.step=20301-100;
[m,prior]=sippi_prior(prior);
%m_old=m;
%[m,prior]=sippi_prior(prior,m);
% imagesc(m{1}-m_old{1});colorbar;axis image;caxis([-1 1].*1e-5)
figure(13);clf
for i=1:5;
    m_old=m;
    [m,prior]=sippi_prior(prior,m);
    subplot(2,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    ax=axis;
    subplot(2,5,5+i);
    d_diff=m{1}.*0;
    ii=find(m_old{1}==m{1});
    plot(prior{1}.xx(ii),prior{1}.yy(ii),'k.');
    set(gca,'ydir','revers')
    axis(ax)
    
    
    %xlabel('X')
    %ylabel('Y')
end
try;load cmap;colormap(cmap);end
subplot(2,5,5);colorbar_shift;
print_mul('prior_reals_visim_seqgibbs_type2');



