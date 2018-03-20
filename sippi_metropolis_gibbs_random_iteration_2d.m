% sippi_metropolis_gibbs_random_iteration(C,mcmc,i);
%    computes a 2D marginal (using Nm ranbdom voronois cells)
%    and movies to a new realizations
%
% See also: sippi_metropolis_gibbs_random_iteration
%
function [C,mcmc]=sippi_metropolis_gibbs_random_iteration_2d(C,mcmc,i);


mcmc.gibbs.null='';
if ~isfield(mcmc.gibbs,'Nm')
    mcmc.gibbs.Nm=200;
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
        
        m_arr=zeros(2,mcmc.gibbs.Nm);
        logL=zeros(1,mcmc.gibbs.Nm);
        pPrior=zeros(1,mcmc.gibbs.Nm);
        for im=1:mcmc.gibbs.Nm;
            [m_test{im},prior_test{im}]=sippi_prior(prior,m);
            [d_test{im},forward]=sippi_forward(m_test{im},C{ic}.forward,prior_test{im},C{ic}.data);
            [logL(im)]=sippi_likelihood(d_test{im},C{ic}.data);
            
            %m_test{im};
            m_arr(:,im) = [m_test{im}{ip}]';
            
        end
        %compute prior
        for iip=1:length(ip);
            if (strcmp(lower(C{ic}.prior_current{ip(iip)}.type),'uniform'))
                logPrior(iip,:)=ones(1,mcmc.gibbs.Nm)./mcmc.gibbs.Nm;
            elseif (strcmp(lower(C{ic}.prior_current{ip(iip)}.type),'gaussian'))
                logPrior(iip,:) =  log(normpdf(m_arr(iip,:),C{ic}.prior_current{ip(iip)}.m0,C{ic}.prior_current{ip(iip)}.std));
            end
        end
        logPost = logL + sum(logPrior);
        pPost=exp(logPost-max(logPost));
        
        %% Nearest Neighbor interpolate
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
        P=exp(i_logPost_G-max(i_logPost_G(:)));
        
        
        % sample from 2D pdf
        [m1_sim,m2_sim,ir_1,ir_2]=sample_from_2d_pdf(P,m1,m2);
        % find the corresponding prior real
        im_use = find(logPost==i_logPost_G(ir_2,ir_1));
        
        
        
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
        
        %% plot
        if (i/mcmc.i_plot)==round(i/mcmc.i_plot)
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
            plot(m1(ir_1),m2(ir_2),'r*')
            hold off
            xlabel(sprintf('m_%d',ip(1)))
            ylabel(sprintf('m_%d',ip(2)))
            title(sprintf('f(m_%d,m_%d | not) at iteration %d',ip(1),ip(2),mcmc.i))
            
            
            drawnow;
        end
    else
        sippi_verbose(sprintf('%s: cannot perfomr Gibbs sampling on selected priors',mfilename));
    end
    
    
end