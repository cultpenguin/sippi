clear all; close all

use_reference=0;
use_prior=6; 
use_metropolis=1;
use_rejection=0;
n_ite=400000;
n_ite=400000;
n_reals_out=500;
doAnneal=1; 
doTempering=0;
dx=0.2;

%plot_posterior_sample=0;
%plot_posterior_data=0;
%plot_posterior_loglikelihood=1;
plot_posterior=0;

n_runs=3;
for i=1:n_runs;
    
    clear prior* options*;
    rseed=i+1;
    sippi_AM13;
    o{i}=options;
    
end
%rseed=11;
%sippi_AM13;
%o{2}=options;%

txt=sprintf('SIPPI_A13_test_mixing_P%d_REF%d',use_prior,use_reference);
save(txt);


%% 
%load SIPPI_A13_test_mixing_P2;
%load SIPPI_A13_test_mixing_P4_REF1.mat
skip_seq_gibbs=0;
for io=1:length(o);
    [reals{io},etype_mean{io},etype_var{io},reals_all{io},reals_ite{io}]=sippi_get_sample(o{io}.txt,1,15,skip_seq_gibbs);
end
close all 

%% IMAGE ETYPE + REF
figure(2);clf;
for io=1:length(o);
    subplot(1,length(o)+1,io+1);
    imagesc(prior{1}.x,prior{1}.y,etype_mean{io});
    axis image
    caxis(prior{1}.cax)    
    title(sprintf('RUN #%02d',io))
end
if isfield(options.mcmc,'m_ref')
    subplot(1,length(o)+1,1);
    imagesc(prior{1}.x,prior{1}.y,options.mcmc.m_ref{1});
    caxis(prior{1}.cax)
    title('Ref')
end
axis image
print_mul(sprintf('%s_etype',txt))   

%% COMPUTE CC
nr=size(reals_all{1},1);
i_sample=o{1}.mcmc.i_sample;
nite=o{1}.mcmc.nite;
i_ax=i_sample:i_sample:nite;


col{1}=[0 0 0];
col{2}=[1 0 0];
col{3}=[0 0 1];
lt{1}='-';
lt{2}=':';
lt{3}='--';


clear h* cc* leg*
ih=0;
figure(3);clf;
for io=1:length(o);;
    
    for i=1:nr
        if isfield(options.mcmc,'m_ref')
            cc=corrcoef(options.mcmc.m_ref{1}(:),reals_all{io}(i,:)');
            cc_ref{io}(i)=cc(2);
        end
        
        % i_use=ceil(nr/3);
        i_use=nr;
        for io_end=1:length(o);;
            cc=corrcoef(reals_all{io_end}(i_use,:),reals_all{io}(i,:)');
            cc_mix{io}{io_end}(i)=cc(2);
        end
    end
    
    if isfield(options.mcmc,'m_ref')
        plot(i_ax,cc_ref{io},'-','color',col{io},'LineWidth',2);
        ih=ih+1;leg{ih}=sprintf('Ref_%d',io)
        drawnow,
        if ih==1; hold on; end
    end
    for io_end=1:length(o);;
        if io_end==io; 
            lw=3;
        else
            lw=1;
        end
        plot(i_ax,cc_mix{io}{io_end},'-','color',col{io},'LineWidth',lw);
        ih=ih+1;leg{ih}=sprintf('mix_{%d,%d}',io,io_end)
        if ih==1; hold on; end
    end
    
    
        
end
hold off
legend(leg,'Location','NorthEastOutside')
print_mul(sprintf('%s_mixing',txt))   


