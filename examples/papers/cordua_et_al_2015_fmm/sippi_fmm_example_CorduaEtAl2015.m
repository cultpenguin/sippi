% sippi_fmm_example_CorduaEtAl2015: CorduaEtAl2015 Figure 4, and Figure 6
%
% Example replicating the example in section 3 in 
%  "
%  Cordua et al., 2015 - Improving the Pattern Reproducibility of
%  Multiple-Point-Based Prior Models Using Frequency Matching.
%  Mathematical Geosciences 47(3).
% "
%
% See also: sippi_forward_fmm, sippi_likelihood_fmm, multinomial
%
clear all;close all
options.txt='CorduaEtAl2015'; % string to append to output files

%% forward 
forward.forward_function='sippi_forward_fmm';
forward.fmm.N=9;

%% prior;
ip=1;
prior{ip}.type='fftma';
dx=0.25;
prior{ip}.x=0:dx:5;
prior{ip}.y=0:dx:12;
prior{ip}.m0=0;

% Correlated prior
%prior{ip}.Cm='1 Gau(3,90,0.25)'; % horizontal
%prior{ip}.d_target=[0 0 0 0 0 0 0 1 1 1];
% uniform prior
prior{ip}.Cm='1 Sph(.5,90,0.2)';
prior{ip}.d_target=[0 1];

[m_ref,prior]=sippi_prior(prior);



%% REFERENCE DATA
use_ti=1;
if use_ti==1;
    TI=read_eas('thinnedti.dat');
    TI=reshape(TI,84,84)';
    %TI=channels;
    %j=2;
    %TI=TI(j:j:end,j:j:end);
    m_ref{1}=TI;
    x_ti=[1:1:size(TI,2)].*dx;
    y_ti=[1:1:size(TI,1)].*dx;
end
[d,forward]=sippi_forward_fmm(m_ref,forward);
data{1}.d_obs=d{1};
data{1}.noise_model='sippi_likelihood_fmm';
data{1}.nprior=5; % counts in the prior uniform fmm distribution

%%
figure(2);clf;
subplot(1,2,1);
imagesc(x_ti,y_ti,TI);axis image;
subplot(1,2,2);
m_ex=sippi_prior(prior);
imagesc(prior{1}.x,prior{1}.y,m_ex{1});
axis image

%% test prior/forward/likelihood
m=sippi_prior(prior);
[d,forward]=sippi_forward(m,forward);
[logL]=sippi_likelihood(d,data);

%% MCMC
prior{1}.seq_gibbs.type=2;
prior{1}.seq_gibbs.i_update_step_max=2000;
prior{1}.seq_gibbs.step=2000;
prior{1}.seq_gibbs.step_min=2;
prior{1}.seq_gibbs.step_max=2000;;

options.mcmc.nite=10000;   % [1] : Number if iterations
options.mcmc.i_sample=ceil(options.mcmc.nite/100); % : Number of iterations between saving model to disk
options.mcmc.i_plot=1000;  % [1]: Number of iterations between updating plots
profile off
t_start=now;
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
t_mcmc=(now-t_start)*3600*24;
save(options.txt);

%% FMM ANALYSIS
clear H*
Nmax = 20; % image stat for the Nmax most frequent TI patterns
sdata=flipud(sort(data{1}.d_obs(:)));
maxd=sdata(min([Nmax,length(data{1}.d_obs)]));
ibin=find(data{1}.d_obs>=maxd);
ibin=ibin(1:Nmax);
H_bin_ti=data{1}.d_obs(ibin)./sum(data{1}.d_obs);

n_reals=5; % number of 'posterior' realisations to consider..
skip_burn_in=1;
m_post=sippi_get_sample(options.txt,1,n_reals,skip_burn_in);
for i=1:size(m_post,3)
    m_c{1}=m_post(:,:,i);    
    [d,forward]=sippi_forward_fmm(m_c,forward);
    H_bin(:,i)=d{1}(ibin)./sum(d{1});
    H_bin_diff(:,i)=H_bin_ti-H_bin(:,i);
    if i<=5
        figure(11);
        subplot(1,5,i);
        imagesc(prior{1}.x,prior{1}.y,m_c{1});;
        axis image
    end
end
figure(11);set_paper;
print_mul(sprintf('%s_posterior_reals',options.txt));

figure(4);clf;set_paper;
subplot(1,2,1)
bar(1:Nmax,[H_bin_ti,H_bin]);
legend('TI','real #1');
xlabel(sprintf('%d most probable TI patterns',Nmax));
ylabel('Frequency')

subplot(1,2,2)
bar(1:Nmax,[H_bin_diff]);
legend('real #1','real #2');
xlabel(sprintf('%d most probable TI patterns',Nmax));
ylabel('Difference to TI freq')

print_mul(sprintf('%s_fmm_match',options.txt));







