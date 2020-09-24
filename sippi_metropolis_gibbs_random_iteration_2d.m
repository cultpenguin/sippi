% sippi_metropolis_gibbs_random_iteration(C,mcmc,i);
%    computes a 2D marginal (using Nm ranbdom voronois cells)
%    and movies to a new realizations
%
% See also: sippi_metropolis_gibbs_random_iteration
%
function [C,mcmc]=sippi_metropolis_gibbs_random_iteration_2d(C,mcmc,i,doPlot);


mcmc.gibbs.null='';
if ~isfield(mcmc.gibbs,'Nm')
    mcmc.gibbs.Nm=200;
end

if ~isfield(mcmc.gibbs,'useNN')
    mcmc.gibbs.useNN=1;
end
if nargin<3, i=1; end
if nargin<4, 
    if (i/mcmc.i_plot)==round(i/mcmc.i_plot)
        doPlot=1;
    else 
        doPlot=0;
    end
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
    
    
    if length(usep)>1
        
        % choose a random set of 2D marginal
        ip=randomsample(usep,2);
        sippi_verbose(sprintf('%s: Running Gibbs sampling of im=[%d,%d] at iteration %d',mfilename,ip(1),ip(2),mcmc.i),0)
        
        prior = C{ic}.prior_current;
        m = C{ic}.m_current;
        
        % update which prior parameters to use
        for j=1:length(prior);
            if sum(ip==j)
                prior{j}.perturb=1;
                prior{j}.seq_gibbs.step=1;
            else
                prior{j}.perturb=0;
                prior{j}.seq_gibbs.step=1;
            end
        end
        
        
        % generate mcmc.gibbs.Nm realizations of the prior, and compute
        % the corresponding likelihood of each.
        m_arr=zeros(2,mcmc.gibbs.Nm);
        logL=zeros(1,mcmc.gibbs.Nm);
        pPrior=zeros(1,mcmc.gibbs.Nm);
        for im=1:mcmc.gibbs.Nm;
            if im==1;
                % make sure current model is part of the tested models
                m_test{im}=C{ic}.m_current;
                prior_test{im}=C{ic}.prior_current;
            else
                [m_test{im},prior_test{im}]=sippi_prior(prior,m);
            end
            %[m_test{im},prior_test{im}]=sippi_prior(prior,m);
            [d_test{im},forward]=sippi_forward(m_test{im},C{ic}.forward,prior_test{im},C{ic}.data);
            [logL(im)]=sippi_likelihood(d_test{im},C{ic}.data);
            
            %m_test{im};
            m_arr(:,im) = [m_test{im}{ip}]';
            
        end
        
        % compute the prior probability 
        for iip=1:length(ip);
            if (strcmp(lower(C{ic}.prior_current{ip(iip)}.type),'uniform'))
                logPrior(iip,:)=ones(1,mcmc.gibbs.Nm)./mcmc.gibbs.Nm;
            elseif (strcmp(lower(C{ic}.prior_current{ip(iip)}.type),'gaussian'))
                logPrior(iip,:) =  log(normpdf(m_arr(iip,:),C{ic}.prior_current{ip(iip)}.m0,C{ic}.prior_current{ip(iip)}.std));
            end
        end
        
        if isfield(mcmc,'anneal')
            [C{ic}.T_fac,mcmc]=sippi_anneal_temperature(i,mcmc,C{ic}.prior_current);
            T=C{ic}.T_fac.*C{ic}.T;
        else
            T=C{ic}.T; % only use the 'base' temperature.
        end
        logPost = logL + sum(logPrior);
        logPost = logL;
        logPost_norm = logPost-max(logPost);
        
        
        if mcmc.gibbs.useNN==1
            %% Nearest Neighbor interpolate
            % 1) Interpolate to 2D conditional probabiility
            % 2) generate a relazation from the 2D
            % 3) find the the precompute set of [m1,m2] closest to the
            % simulated set an dchoose that as a realizations!
            % Sometimes this measn the chain will move to a place with poor
            % data fit!
            
            if isfield(mcmc.gibbs,'Nn1');
                Nn1 = mcmc.gibbs.Nn1;
            else
                Nn1=31;
            end
            if isfield(mcmc.gibbs,'Nn2');
                Nn2 = mcmc.gibbs.Nn2;
            else
                Nn2=30;
            end
            m1 = linspace(min(m_arr(1,:)),max(m_arr(1,:)),Nn1);
            m2 = linspace(min(m_arr(2,:)),max(m_arr(2,:)),Nn2);
            im1=interp1([m1(1) m1(end)],[1 Nn1],m_arr(1,:));
            im2=interp1([m2(1) m2(end)],[1 Nn2],m_arr(2,:));
            
            [imm1,imm2]=meshgrid(1:Nn1,1:Nn2);
            i_logPost_G=griddata(im1,im2,logPost,imm1,imm2,'nearest');
            %logPost_G=griddata(m_arr(1,:),m_arr(2,:),logPost,mm1,mm2,'nearest');
            P=exp( (1/T)* (i_logPost_G-max(i_logPost_G(:))) );
            
            %keyboard
            % sample from 2D pdf
            [m1_sim,m2_sim,ir_1,ir_2]=sample_from_2d_pdf(P,m1,m2);
            
            % Find closest match
            w=max(m_arr')-min(m_arr');
            %dx = max
            for im=1:mcmc.gibbs.Nm;
                dis(im) = sqrt( ((m_test{im}{ip(1)}-m1_sim)./(w(1)))^2 + ((m_test{im}{ip(2)}-m2_sim)./(w(2)))^2 );
            end
            im_use=find(dis==min(dis));
            im_use=im_use(1);
        else
            %% USE A RANDOM ONE FO THE PROPOSED MODELS
            P_acc = exp((1/T)*(logPost-max(logPost)));
            Cum_P = cumsum(P_acc);
            Cum_P = Cum_P./max(Cum_P);
            dp=1/mcmc.gibbs.Nm;
            p=[1:mcmc.gibbs.Nm]*dp;
            
            r=rand(1);
            %r=0.04
            im_use = min(find(Cum_P>r));
            im_use=im_use(1);
            
        end
        
        
        % Next line only for debugging        
        logL_old = C{ic}.logL_current;
        
        
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
        C{ic}.L_current=logL(im_use);
        C{ic}.iacc=C{ic}.iacc+1;
        %C{ic}.mcmc.logL(C{ic}.iacc)=C{ic}.logL_current;
        C{ic}.mcmc.acc(ip,mcmc.i)=1;
        C{ic}.mcmc.logL(mcmc.i)=C{ic}.logL_current;
        
        
        % Next 6 lines only for debugging        
        logL_new = C{ic}.logL_current;
        %P_acc = min([2 exp(logL_new-logL_old)]);
        %disp(sprintf('logL New/Old  %5.3f/%5.3f, Pacc=%1.11f',logL_new, logL_old,P_acc))
        %if P_acc<0.00001
        %    doPlot=1;
        %    keyboard
        %end
        
        
        
        %% plot
        if mcmc.gibbs.useNN==1
            if (doPlot==1)
                figure_focus(63);
                subplot(1,2,1);
                imagesc(m1,m2,i_logPost_G);
                hold on
                plot(m_arr(1,:),m_arr(2,:),'r.');
                plot(m1(ir_1),m2(ir_2),'g*')
                hold off
                xlabel(sprintf('m_%d',ip(1)))
                ylabel(sprintf('m_%d',ip(2)))
                title(sprintf('log(f(m_%d,m_%d | not)) at iteration %d',ip(1),ip(2),mcmc.i))
                
                subplot(1,2,2);
                imagesc(m1,m2,P);colormap(1-gray)
                hold on
                plot(m1(ir_1),m2(ir_2),'r.')
                plot([m{ip(1)},C{ic}.m_current{ip(1)}],[m{ip(2)},C{ic}.m_current{ip(2)}],'r--')
                plot(C{ic}.m_current{ip(1)},C{ic}.m_current{ip(2)},'r*')
                plot(m{ip(1)},m{ip(2)},'g*')
                hold off
                xlabel(sprintf('m_%d',ip(1)))
                ylabel(sprintf('m_%d',ip(2)))
                title(sprintf('f(m_%d,m_%d | not) at iteration %d, T=%g',ip(1),ip(2),mcmc.i, T))
                
                drawnow;
            end
        else
            if (doPlot==1)
                %%
                figure_focus(63);
                subplot(1,2,1);
                plot(m_arr(1,:),m_arr(2,:),'k.','MarkerSize',10);hold on;;scatter(m_arr(1,:),m_arr(2,:),12,logPost,'filled');hold off
                hold on
                plot(m_arr(1,im_use),m_arr(2,im_use),'b*');
                plot(m{ip(1)},m{ip(2)},'g*');
                hold off
                xlabel(sprintf('m_%d',ip(1)))
                ylabel(sprintf('m_%d',ip(2)))
                title(sprintf('log(f(m_%d,m_%d | not)) at iteration %d',ip(1),ip(2),mcmc.i))
                colormap(gca,flipud(hot));
                colorbar
                
                subplot(1,2,2);
                plot(m_arr(1,:),m_arr(2,:),'.','MarkerSize',1,'color',[.5 .5 .5]);
                hold on;;
                try
                scatter(m_arr(1,:),m_arr(2,:),P_acc*100+.01,P_acc,'filled');
                end
                plot(m_arr(1,im_use),m_arr(2,im_use),'b*');
                plot(m{ip(1)},m{ip(2)},'g*');
                hold off
                colormap(gca,flipud(hot));
                colorbar
                xlabel(sprintf('m_%d',ip(1)))
                ylabel(sprintf('m_%d',ip(2)))
                title(sprintf('f(m_%d,m_%d | not) at iteration %d, T=%g',ip(1),ip(2),mcmc.i, T))
                %hold on
                drawnow;
            end
        end
        
    else
        sippi_verbose(sprintf('%s: cannot perfomr Gibbs sampling on selected priors',mfilename));
    end
    
    
end