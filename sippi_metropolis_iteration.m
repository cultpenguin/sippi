function [C,mcmc]=sippi_metropolis_iteration(C,mcmc,i)

NC=length(C);
for ic=1:NC
    
    %% set seed / necessary?
    for im=1:length(C{ic}.prior_current)
        C{ic}.prior_current{im}.seed=i;
    end
    
    %% SELECT IF STEP LENGTH HAS TO BE UPDATED
    for im=1:length(C{ic}.prior_current)
        if (( (mcmc.i./C{ic}.prior_current{im}.seq_gibbs.i_update_step)==round((mcmc.i./C{ic}.prior_current{im}.seq_gibbs.i_update_step)))&(mcmc.i<C{ic}.prior_current{im}.seq_gibbs.i_update_step_max))
            % UPDATE STEP LENGTH
            C{ic}.prior_current=sippi_prior_set_steplength(C{ic}.prior_current,C{ic}.mcmc,im);
        end
        C{ic}.mcmc.step(im,i)=C{ic}.prior_current{im}.seq_gibbs.step(1);
    end
    
    %% Sample prior
    % SELECT WHICH MODEL PARAMETERS TO PERTURB    
    for im=1:length(C{ic}.prior_current);
        C{ic}.prior_current{im}.perturb=0;
    end
    
    % PERTURBATION STRATEGY
    if mcmc.pert_strategy.perturb_all==1,
        % perturb all
        im_perturb=1:1:length(C{ic}.prior_current);
    elseif mcmc.pert_strategy.perturb_all==2,
        % perturb random selection
        i_perturb=round(rand(1,length(C{ic}.prior_current)));
        im_perturb=find(i_perturb);
        if isempty(im_perturb);
            im_perturb=ceil(rand(length(C{ic}.prior_current)));
        end
    else
        % perturb one parameter according to frequency distribution
        i_pert=mcmc.pert_strategy.i_pert;
        pert_freq=cumsum(mcmc.pert_strategy.i_pert_freq);
        pert_freq=pert_freq./max(pert_freq);
        im_perturb=i_pert(min(find(rand(1)<pert_freq)));       
    end
    for k=1:length(im_perturb);
        C{ic}.prior_current{im_perturb(k)}.perturb=1;
        C{ic}.mcmc.perturb(im_perturb(k),i)=1;
    end
    
    % SAMPLE PRIOR
    [C{ic}.m_propose,C{ic}.prior_propose] = sippi_prior(C{ic}.prior_current,C{ic}.m_current);
    
    %% FORWARD PROBLEM
    [C{ic}.d,C{ic}.forward,C{ic}.prior_propose,C{ic}.data]=sippi_forward(C{ic}.m_propose,C{ic}.forward,C{ic}.prior_propose,C{ic}.data);
    
    %% LIKELIHOOD
    [C{ic}.logL_propose,C{ic}.L_propose,C{ic}.data]=sippi_likelihood(C{ic}.d,C{ic}.data);
    
    %%  MOVE?
    % Accept probability
    
    % set temperature
    if mcmc.do_anneal==1;
        [C{ic}.T_fac,mcmc]=sippi_anneal_temperature(i,mcmc,C{ic}.prior_current);
        T(ic)=C{ic}.T_fac.*C{ic}.T;
    else
        T(ic)=C{ic}.T;
    end
    
    C{ic}.Pacc = exp((1./T(ic)).*(C{ic}.logL_propose-C{ic}.logL_current));
    
    if (mcmc.accept_only_improvements==1)
        % Optimization only?
        if C{ic}.logL_propose>C{ic}.logL_current
            C{ic}.Pacc=1;
        else
            C{ic}.Pacc=0;
        end
    end
    % Optionally accept all proposed models
    if (mcmc.accept_all==1), C{ic}.Pacc=1; end
    
    % Accept move
    C{ic}.forward.last_proposed_model_accept=0; %
    if C{ic}.Pacc>rand(1)
        % ACCEPT MODEL
        C{ic}.forward.last_proposed_model_accept=1;
        
        % move to current model
        C{ic}.prior_current=C{ic}.prior_propose; % NEEDED FOR GAUSSIAN TYPE PRIOR
        C{ic}.m_current=C{ic}.m_propose;
        C{ic}.d_current=C{ic}.d;
        C{ic}.logL_current=C{ic}.logL_propose;
        C{ic}.L_current=C{ic}.L_propose;
        C{ic}.iacc=C{ic}.iacc+1;
        %C{ic}.mcmc.logL(C{ic}.iacc)=C{ic}.logL_current;
        C{ic}.mcmc.acc(im_perturb,mcmc.i)=1;
        
    else
        % REJECT MOVE MODEL
    end
    C{ic}.mcmc.logL(mcmc.i)=C{ic}.logL_current;
    
    
    if mcmc.store_all==1
        %C{ic}.logL_all(i)=C{ic}.logL_current;
        for im=1:length(prior);
            C{ic}.m_cur_all(im,i) = C{ic}.m_current{im}(1);
        end
    end
    
end % END MC CHAIN LOOP

%% TEMPERING/ORDERING
if length(C)>1
    if rand(1)<mcmc.chain_frequency_jump; %frequency of swap tests
        % IF TEMPERING
        ic_i=ceil(NC*rand(1));
        j_arr=setxor(1:NC,ic_i);
        ic_j=j_arr(ceil((NC-1)*rand(1)));
        
        Pswap = exp((C{ic_j}.logL_current)*(1./T(ic_i) - 1./T(ic_j)) + ...
                    (C{ic_i}.logL_current)*(1./T(ic_j) - 1./T(ic_i)) );
        
        if rand(1)<Pswap
            % accept swap
            
            mcmc.n_swap=mcmc.n_swap+1;
            
            mcmc.i_swap(1,mcmc.n_swap)=ic_i;
            mcmc.i_swap(2,mcmc.n_swap)=ic_j;
            
            % perform the swap
            C_i=C{ic_i};
            
            C{ic_i}.m_current=C{ic_j}.m_current;
            C{ic_i}.prior_current=C{ic_j}.prior_current;
            C{ic_i}.data=C{ic_j}.data;
            C{ic_i}.logL_current=C{ic_j}.logL_current;
            C{ic_i}.m_current=C{ic_j}.m_current;
            C{ic_i}.i_chain=C{ic_j}.i_chain;
            
            C{ic_j}.m_current=C_i.m_current;
            C{ic_j}.prior_current=C_i.prior_current;
            C{ic_j}.data=C_i.data;
            C{ic_j}.logL_current=C_i.logL_current;
            C{ic_j}.m_current=C_i.m_current;
            C{ic_j}.i_chain=C_i.i_chain;
            
            % Keep step length constant within chains
            for k=1:NC;for im=1:length(C{k}.prior_current);
                    C{k}.prior_current{im}.seq_gibbs.step=C{k}.mcmc.step(im,i);
                end;end
            sippi_verbose(sprintf('%s: at i=%05d SWAP chains [%d<->%d]',mfilename,i,ic_i,ic_j),2);
        end
    end
    
    if (mcmc.order_chains_frequency>0)
        if (rand(1)<mcmc.order_chains_frequency);
            for ic=1:NC
                logL(ic)=C{ic}.logL_current;
            end
            
            % select best chain from likelihood
            P=exp(logL-max(logL));
            cP=cumsum(P)/sum(P);
            i1=find(cP>rand(1));
            ic_1=i1(1);
            if ic_1~=1
                sippi_verbose(sprintf('Ordering chain. Using i_chain=%d',C{ic_1}.i_chain),1)
                ic_i=1;
                ic_j=ic_1;
                
                % perform the swap
                C_i=C{ic_i};
                
                C{ic_i}.m_current=C{ic_j}.m_current;
                C{ic_i}.prior_current=C{ic_j}.prior_current;
                C{ic_i}.data=C{ic_j}.data;
                C{ic_i}.logL_current=C{ic_j}.logL_current;
                C{ic_i}.m_current=C{ic_j}.m_current;
                C{ic_i}.i_chain=C{ic_j}.i_chain;
                
                C{ic_j}.m_current=C_i.m_current;
                C{ic_j}.prior_current=C_i.prior_current;
                C{ic_j}.data=C_i.data;
                C{ic_j}.logL_current=C_i.logL_current;
                C{ic_j}.i_chain=C_i.i_chain;
                
            end
        end
    end
    
    %% store chain id
    for ic=1:NC
        mcmc.i_chain(ic,i)=C{ic}.i_chain;
    end
    
end

