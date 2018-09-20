function [C,mcmc]=sippi_metropolis_gibbs_random_iteration(C,mcmc,i);


mcmc.gibbs.null='';
if ~isfield(mcmc.gibbs,'N_bins')
    mcmc.gibbs.N_bins=31;
end
if isfield(mcmc.gibbs,'Nm')
    mcmc.gibbs.N_bins=mcmc.gibbs.Nm;
end

NC=length(C);
for ic=1:NC;
    
    % CHECK NUMBER OF 1D PRIORS
    if isfield(mcmc.gibbs,'i_pert');
        usep = mcmc.gibbs.i_pert;
    else
        usep=[];
        for ip=1:length(C{ic}.prior_current)
            if (C{ic}.prior_current{ip}.ndim==0)
                usep=[usep ip];
            end
        end
    end
    
    if length(usep)>0
        ip=ceil(rand(1)*length(usep));
        sippi_verbose(sprintf('%s: Running Gibbs sampling of im=%d at iteration %d',mfilename,ip,mcmc.i),1)
        
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
        
        m_arr=zeros(1,mcmc.gibbs.N_bins);
        logL=zeros(1,mcmc.gibbs.N_bins);
        pPrior=zeros(1,mcmc.gibbs.N_bins);
        for im=1:mcmc.gibbs.N_bins;
            if im==1;
                % make sure current model is part of the tested models
                m_test{im}=C{ic}.m_current;
                prior_test{im}=C{ic}.prior_current;
            else
                [m_test{im},prior_test{im}]=sippi_prior(prior,m);
            end
            [d_test{im},forward]=sippi_forward(m_test{im},C{ic}.forward,prior_test{im},C{ic}.data);
            [logL(im)]=sippi_likelihood(d_test{im},C{ic}.data);
            
            %m_test{im};
            m_arr(im) = m_test{im}{ip};
            
        end
        %compute prior
        if (strcmp(lower(C{ic}.prior_current{ip}.type),'uniform'))
            logPrior=ones(1,mcmc.gibbs.N_bins)./mcmc.gibbs.N_bins;
        elseif (strcmp(lower(C{ic}.prior_current{ip}.type),'gaussian'))
            logPrior =  log(normpdf(m_arr,C{ic}.prior_current{ip}.m0,C{ic}.prior_current{ip}.std));
        end
        
        % get the posterio marginal
        logPost = logL + logPrior;
        % anealing type probability
        % if mcmc.do_anneal==1 % faster, bust introdcuces an extra variable
        if isfield(mcmc,'anneal')
            [C{ic}.T_fac,mcmc]=sippi_anneal_temperature(i,mcmc,C{ic}.prior_current);
            T=C{ic}.T_fac.*C{ic}.T;
        else
            T=C{ic}.T; % only use the 'base' temperature.
        end
        %pPost=exp(logPost-max(logPost)).^(1/T);
        pPost=exp( (1/T)*(logPost-max(logPost)) );
        
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
        im_use=im_use(1); % just in case the prior generates the same realization
        
        
        % move to current model
        prior_ref=C{ic}.prior_current;
        C{ic}.prior_current=prior_test{im_use}; % NEEDED FOR GAUSSIAN TYPE PRIOR
        % keep step length (all pars in seq Gibbs)
        for im=1:length(prior_ref);
            C{ic}.prior_current{im}.seq_gibbs=prior_ref{im}.seq_gibbs;
        end
        
        C{ic}.m_current=m_test{im_use};
        C{ic}.d_current=d_test{im_use};
        C{ic}.logL_current=logL(im_use);
        C{ic}.logL_propose=logL(im_use);
        C{ic}.L_current=logL(im_use);
        C{ic}.iacc=C{ic}.iacc+1;
        %C{ic}.mcmc.logL(C{ic}.iacc)=C{ic}.logL_current;
        C{ic}.mcmc.acc(ip,mcmc.i)=1;
        try
            C{ic}.mcmc.logL(mcmc.i)=C{ic}.logL_current;
        catch
            keyboard
        end
        if (i/mcmc.i_plot)==round(i/mcmc.i_plot)
            %% plot
            figure_focus(63);
            subplot(1,3,1);
            plot(m_arr,[logL;logPrior;logPost],'.')
            xlabel(sprintf('m_%d -%s',ip,prior{ip}.name))
            ylabel(sprintf('f(m_%d | f_{not %d})',ip,ip))
            subplot(1,3,2);
            plot(m_arr,pPost,'.',s_m_arr,s_pPost,'r-');
            ylim([0 1])
            hold on
            plot([1 1].*m_sim,ylim,'r-')
            plot([1 1].*m{ip},ylim,'g--')
            hold off
            xlabel(prior{ip}.name)
            title(sprintf('f(m_%d | not m_%d), #ite=%d, T=%g',ip,ip,mcmc.i,T))
            
           
            subplot(1,3,3);
            plot(s_m_arr,s_cpdf,'k.','MarkerSize',4)
            hold on
            plot(m_sim,s_cpdf(i_sim),'r.','MarkerSize',30)
            
            %plot(s_m_arr,s_cpdf,'-*');
            grid on;
            plot(xlim,[1 1].*r);
            plot([1 1].*m_sim,ylim,'r--')
            hold off
            drawnow;
        end
    else
        sippi_verbose(sprintf('%s: cannot perfomr Gibbs sampling on selected priors',mfilename));
    end
    
    
    
    
end