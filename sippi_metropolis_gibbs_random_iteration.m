function [C,mcmc]=sippi_metropolis_gibbs_random_iteration(C,mcmc,i);

        
mcmc.gibbs.null='';
if ~isfield(mcmc.gibbs,'N_bins')
    mcmc.N_bins=31;
end

NC=length(C);
for ic=1:NC;
    
    % CHECK NUMBER OF 1D PRIORS
    usep=[];
    for ip=1:length(C{ic}.prior_current)
        if (C{ic}.prior_current{ip}.ndim==0)
            usep=[usep ip];
        end
    end
    
    if length(usep)>0
        ip=ceil(rand(1)*length(usep));
        
        mcmc.N_bins = 151;
        prior = C{ic}.prior_current;
        m = C{ic}.m_current;
        
        % update which prior parameters to use
        for j=1:length(prior);
            if ip==j
                prior{j}.perturb=1;
                prior{j}.seq_gibbs.step=1;
            else
                prior{j}.perturb=0;
                prior{j}.seq_gibbs.step=1;
            end
        end
        
        m_arr=zeros(1,mcmc.N_bins);
        logL=zeros(1,mcmc.N_bins);
        pPrior=zeros(1,mcmc.N_bins);
        for im=1:mcmc.N_bins;
            [m_test{im},prior_test{im}]=sippi_prior(prior,m);            
            [d_test{im},forward]=sippi_forward(m_test{im},C{ic}.forward,prior_test{im},C{ic}.data);
            [logL(im)]=sippi_likelihood(d_test{im},C{ic}.data);
            
            %m_test{im};
            m_arr(im) = m_test{im}{ip};
            
        end
        %compute prior
        if (strcmp(lower(C{ic}.prior_current{ip}.type),'uniform'))
            logPrior=ones(1,mcmc.N_bins)./mcmc.N_bins;
        elseif (strcmp(lower(C{ic}.prior_current{ip}.type),'gaussian'))
            logPrior =  log(normpdf(m_arr,C{ic}.prior_current{ip}.m0,C{ic}.prior_current{ip}.std));
        end

        logPost = logL + logPrior;
        pPost=exp(logPost-max(logPost));
        
        sd=sortrows([m_arr;pPost]',1);
        s_m_arr = sd(:,1);
        s_pPost = sd(:,2);
        s_cpdf = cumsum(s_pPost);
        s_cpdf = s_cpdf./max(s_cpdf);
        
        % simulated
        r = rand(1);
        i_sim = min(find(s_cpdf>r));
        m_sim=s_m_arr(i_sim);
        
        % update and move to the current model
        im_use = find(m_sim==m_arr);
        
        
         % move to current model
        C{ic}.prior_current=prior_test{im_use}; % NEEDED FOR GAUSSIAN TYPE PRIOR
        C{ic}.m_current=m_test{im_use};
        C{ic}.d_current=d_test{im_use};
        C{ic}.logL_current=logL(im_use);        
        C{ic}.L_current=logL(im_use);
        C{ic}.iacc=C{ic}.iacc+1;
        %C{ic}.mcmc.logL(C{ic}.iacc)=C{ic}.logL_current;
        C{ic}.mcmc.acc(ip,mcmc.i)=1;
        C{ic}.mcmc.logL(mcmc.i)=C{ic}.logL_current;
        
        %% plot
        if (i/mcmc.i_plot)==round(i/mcmc.i_plot)
        figure_focus(63);
        subplot(1,3,1);        
        plot(m_arr,[logL;logPrior;logPost],'.')
        xlabel(sprintf('m_%d',ip))
        ylabel(sprintf('f(m_%d | f_{not %d})',ip,ip))
        subplot(1,3,2);
        plot(m_arr,pPost,'.',s_m_arr,s_pPost,'r-')
        subplot(1,3,3);
        plot(s_m_arr,s_cpdf,'k.','MarkerSize',30)
        hold on
        plot(m_sim,s_cpdf(i_sim),'r.','MarkerSize',30)
        
        hold off
        
        drawnow;
        end
    else
            sippi_verbose(sprintf('%s: cannot perfomr Gibbs sampling on selected priors',mfilename));
    end
    
    
    
    
end