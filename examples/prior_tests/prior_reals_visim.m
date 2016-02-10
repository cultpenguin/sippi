% prior_reals_visim Sampling a VISIM type prior model

clear all,close all
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

rng(1);
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
colorbar_shift;
print_mul('prior_reals_visim');



%%
N=10000;
prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
prior{im}.d_target=d_target;
rng(1);
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
colorbar_shift;
print_mul('prior_reals_visim_target');

%% SEQ GIBBS TYPE 1
try;prior{im}=rmfield(prior{im},'d_target');end
rng(2)
prior{im}.seq_gibbs.type=1;
prior{im}.seq_gibbs.step=[4];

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
colormap(sippi_colormap(1));
subplot(2,5,5);colorbar_shift;
print_mul('prior_reals_visim_seqgibbs_type1');


%% SEQ GIBBS TYPE 2
try;prior{im}=rmfield(prior{im},'d_target');end
rng(1);
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
colormap(sippi_colormap(1));
subplot(2,5,5);colorbar_shift;
print_mul('prior_reals_visim_seqgibbs_type2');


%% Conditional simulation SGSIM type
figure(14);
prior{1}.method='dssim';
prior{1}.d_target=d_target;
[d_pdf,d_x]=hist(d_target,40);
ddx=(d_x(2)-d_x(1));
d_pdf=d_pdf/(sum(d_pdf)*ddx);;
mm_dssim=[];
for i=1:5;
    subplot(2,6,i);
    [m,prior]=sippi_prior(prior);
    mm_dssim=[mm_dssim;m{1}(:)];
    imagesc(prior{1}.x,prior{1}.y,m{1});;axis image
    caxis([8 12])
end
subplot(2,6,6);
[d_pdf_dssim]=hist(mm_dssim,d_x);
d_pdf_dssim=sum(d_pdf)*d_pdf_dssim./sum(d_pdf_dssim);
bar(d_x,d_pdf_dssim);
hold on
plot(d_x,d_pdf,'k-')
hold off
title('method=''dssim''')

prior{1}.method='sgsim';
mm_sgsim=[];
for i=1:5;
    subplot(2,6,6+i);
    [m,prior]=sippi_prior(prior);
    mm_sgsim=[mm_sgsim;m{1}(:)];
    imagesc(prior{1}.x,prior{1}.y,m{1});axis image;
    caxis([8 12])
end
subplot(2,6,12);
[d_pdf_sgsim]=hist(mm_sgsim,d_x);
d_pdf_sgsim=sum(d_pdf)*d_pdf_sgsim./sum(d_pdf_sgsim);
bar(d_x,d_pdf_sgsim);
hold on
plot(d_x,d_pdf,'k-')
hold off
title('method=''sgsim''')
print_mul('prior_reals_visim_sgsim_dssim');


