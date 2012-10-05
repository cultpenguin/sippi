% prior_reals_sisim Sampling a SISIM type prior model
clear all,
im=1; 
prior{im}.type='SISIM';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
prior{im}.Va='1 Sph(10,90,.25)';
prior{im}.marginal_prob=[0.1 0.2 .4 .3];


prior=sippi_prior_init(prior);

randn('seed',1);
%%
figure(10);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
try;load cmap;colormap(cmap);end
colorbar_shift;
print_mul('prior_reals_sisim')
suptitle(sprintf('A priori realizations from a SISIM type prior model'))
