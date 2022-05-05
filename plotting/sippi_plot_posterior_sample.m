% sippi_plot_posterior_sample: plots posterior sample statistics
%
% Call
%    [options]=sippi_plot_posterior_sample(options,prior,data,forward);
%
% See also: sippi_plot_posterior
%
function [options]=sippi_plot_posterior_sample(options,prior,data,forward);

%% LOAD THE CORRECT DATA
cwd=pwd;
if nargin==0
    % LOAD FROM MAT FILES
    [p,matfile]=fileparts(pwd);
    load(matfile);
elseif nargin==1;
    if isstruct(options),
    else
        fname=options;
        cd(fname);
        load(fname);
    end
end

%try
    
    % JUST IN CASE LSQ WAS PERFORMED
    if exist('m_est','var');options.m_est=m_est;end
    if exist('Cm_est','var');;options.Cm_est=Cm_est;end
    
    if nargin<5
        try
            fname=options.txt;
        catch
            fname=mfilename;
        end
    end
    if ~exist('mcmc','var')&isfield(options,'mcmc')
        mcmc=options.mcmc;
        mcmc=options.C{1}.mcmc;
    end
    
    % ALL DATA LOADED
    
    %%
    % SET DFAULT PLOTTING SETTINGS
    options=sippi_plot_defaults(options);
    
    
    %% REALS
    nm=length(prior);
    
    if ~exist('im_arr','var');
        im_arr=1:length(prior);
    end
    
    if ~exist('n_reals','var');
        for j=1:length(im_arr);
            if prior{im_arr(j)}.ndim<2
                n_reals(j)=10000;
            else
                n_reals(j)=5;
            end
        end
    end
    
    if length(n_reals)==1;
        n_reals=ones(1,length(prior)).*n_reals;
    end
    %%
    
    for im=im_arr;
        
        if isfield(prior{im},'name');
            title_txt=sprintf('m%d: %s',im,prior{im}.name);
        else
            title_txt=sprintf('m%d',im);
        end
        
        
        try;cd(plotdir);end
        clear cax;
        % find dimension
        ndim=sum(prior{im}.dim>1);
        %if ndim==0;ndim=1;end
        
        
        % FIND SCALE/ORIENTATION
        ax_lscape=1;
        try
            % if 'lim' is set
            if prior{im}.lim(1)<max(prior{im}.lim(2:3))
                ax_lscape=0;
            end
        end
        try
            % if 'daspect' is set
            r=prior{im}.lim./prior{im}.daspect;
            if r(1)<max(r(2:3))
                ax_lscape=0;
            end
        end
        
        
        options.null='';
        id=1;
        
        
        skip_seq_gibbs=options.plot.skip_seq_gibbs;
        [reals,etype_mean,etype_var,reals_all,ite_reals]=sippi_get_sample('.',im,n_reals(im),skip_seq_gibbs);
        
        m_post{im}=reals;
        
        if ~exist('cax','var');
            if isfield(prior{im},'cax')
                cax=prior{im}.cax;
            else
                try
                    cax=[prior{im}.min prior{im}.max];
                end
            end
        end
        
        x=prior{im}.x;
        y=prior{im}.y;
        z=prior{im}.z;
        
        %% PLOT LAST ACCEPTED MODEL
        try
            if prior{im}.ndim>0
                try
                    sippi_plot_prior(prior,m_current,im);
                    print_mul(sprintf('%s_m%d_last_accepted_model',fname,im),options.plot.hardcopy_types)
                    
                catch
                    sippi_plot_prior(prior,options.mcmc.m_current,im);
                    print_mul(sprintf('%s_m%d_last_accepted_model',fname,im),options.plot.hardcopy_types)
                end
            end
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot last accepted model',mfilename));
            cd(cwd);
        end
        
        
        %% PLOT REFERENCE MODEL
        try
            fn=101;
            figure_focus(fn);clf;
            if isfield(options.mcmc,'m_ref');
                if prior{im}.ndim>0
                    sippi_plot_prior(prior,options.mcmc.m_ref,im,1,fn);
                    watermark(sprintf('reference model - %d - %s',im,prior{im}.name));
                    print_mul(sprintf('%s_m%d_ref_model',fname,im),options.plot.hardcopy_types);
                end
            end
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot reference model',mfilename));
            cd(cwd);
        end
        
        %% PLOT POSTERIOR REALS
        f_id=(im)*10+1;
        figure_focus(f_id);
        set_paper('landscape');clf;
        
        i1=ceil(size(reals,1)/n_reals(im));
        ii=ceil(linspace(i1,size(reals,1),n_reals(im)));
        if ndim==0
            
            %% Plot the value of the 1D prior as a function of iterations number
            try
                figure_focus(im*10+3);
                set_paper('landscape');clf;
                try
                    ii=[1:1:length(reals_all)].*options.mcmc.i_sample;
                catch
                    ii=[1:1:length(reals_all)];
                end
                plot(ii(:),reals_all);
                
                ylim=get(gca,'ylim');
                if isfield(prior{im},'cax');ylim=prior{im}.cax;end
                if isfield(prior{im},'min');ylim(1)=prior{im}.min;end
                if isfield(prior{im},'max');ylim(2)=prior{im}.max;end
                set(gca,'ylim',ylim);
                set(gca,'FontSize',options.plot.axis.fontsize)
                xlabel('Iteration #')
                ylabel(prior{im}.name)
                print_mul(sprintf('%s_m%d_posterior_values',fname,im),options.plot.hardcopy_types)
                
            catch
                disp(sprintf('%s : failed to plot 1D posterior values for prior #%d',mfilename,im));
            end
            
            figure_focus(f_id);
            N=length(reals);
            %N=n_reals(im);
            prior_sample=zeros(1,N);
            clear p;
            p{1}=prior{im};
            for i=1:n_reals(im);
                m=sippi_prior(p);
                sample_prior(i)=m{1};
            end
            
            if ~exist('cax','var');
                cax=[min(sample_prior) max(sample_prior)];
                if isnan(cax(1));cax(1)=0;end
                if isnan(cax(2));cax(2)=0;end
            end
            try
                hx=linspace(cax(1),cax(2),31);
                h_post=hist(reals,hx);
                h_post=h_post/sum(h_post);
                
                h_prior=hist(sample_prior,hx);
                h_prior=h_prior/sum(h_prior);
                
                bar(hx,h_prior,.8,'k');
                hold on
                bar(hx,h_post,.6,'r');
                hold off
                ylim=get(gca,'ylim');
            catch
                sippi_verbose(sprintf('Could not plot bar chart',mfilename),0);
            end
            %% GET p10,050,090
            try
                % ONLY DO IF QUNTILE EXISTS
                p50_post=quantile(reals,.5);
                p_l_post=quantile(reals,.025);
                p_h_post=quantile(reals,.975);
                p50_prior=quantile(sample_prior,.5);
                p_l_prior=quantile(sample_prior,.025);
                p_h_prior=quantile(sample_prior,.975);
                
                hold on
                y0=diff(ylim)*.80+ylim(1);
                yl=diff(ylim)*.74+ylim(1);
                yu=diff(ylim)*.86+ylim(1);
                
                plot([p_l_prior p_h_prior],[1 1].*y0,'k-','LineWidth',3)
                plot([1 1]*p_l_prior,[yl yu],'k-','LineWidth',3)
                plot([1 1]*p_h_prior,[yl yu],'k-','LineWidth',3)
                plot(p50_prior,y0,'k.','MarkerSize',22)
                
                plot([p_l_post p_h_post],[1 1].*y0,'r-','LineWidth',1)
                plot([1 1]*p_l_post,[yl yu],'r-','LineWidth',1)
                plot([1 1]*p_h_post,[yl yu],'r-','LineWidth',1)
                plot(p50_post,y0,'r.','MarkerSize',16)
                hold off
            end
            xlabel(prior{im}.name,'interpreter','none','FontSize',options.plot.axis.fontsize+2)
            %xlabel(prior{im}.name,'FontSize',options.plot.axis.fontsize+2)
            ylabel('Frequency','interpreter','none','FontSize',options.plot.axis.fontsize+2)
            % BUG/20140619 : It seems Matlab R2014b does not handle legend
            % very well when using ppp.m ..
            legend('prior','posterior')
            set(gca,'ytick',[]);
            try
                set(gca,'xlim',cax);
            end
            % PLOT REFERENCE IF IT EXISTS
            try
                if isfield(options.mcmc,'m_ref');
                    hold on
                    %plot(options.mcmc.m_ref{im},y0,'go','MarkerSize',6,'LineWidth',3);
                    plot(options.mcmc.m_ref{im},0,'go','MarkerSize',6,'LineWidth',3);
                    hold off
                end
            end
            
            %ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            
        elseif ndim==1
            
            
            reals = squeeze(reals);
            [nx,nr]=size(reals);
            % check for the used DIM.
            if prior{1}.dim(1)==nx;
                x=prior{1}.x;
            elseif prior{1}.dim(2)==nx;
                x=prior{1}.y;
            elseif prior{1}.dim(2)==nx;
                x=prior{1}.z;
            end
            
            %% SAMPLE DENSITY
            
            figure;
            nbins=21;
            if isfield(prior{im},'cax');
                lim=prior{im}.cax;
                h=linspace(lim(1),lim(2),nbins);
            else
                h=linspace(min(reals(:)),max(reals(:)),nbins);
            end
            nh=length(h);
            
            
            
            probmat=zeros(nh,size(reals,1));
            for i=1:size(reals,1);
                probmat(:,i)=hist(reals(i,:),h);
            end
            probmat=probmat./size(reals,2);
            
            % conditional cdf
            %imagesc(prior{im}.x,h,cumsum(probmat));
            % conditional pdf
            imagesc(x,h,(probmat));
            
            colorbar
            set(gca,'ydir','normal')
            hold on
            colormap(1-gray)
            %colormap(flipud(hot))
            if exist('quantile','file')
                caxis([0 quantile(probmat(:),.90)])
                try
                    plot(x,quantile(squeeze(reals)',.025),'g--','linewidth',1);
                    plot(x,quantile(squeeze(reals)',.5),'g-','linewidth',2);
                    plot(x,quantile(squeeze(reals)',.975),'g--','linewidth',1);
                catch
                    plot(quantile(squeeze(reals)',.025),'g--','linewidth',1);
                    plot(quantile(squeeze(reals)',.5),'g-','linewidth',2);
                    plot(quantile(squeeze(reals)',.975),'g--','linewidth',1);    
                end
            end
            hold off
            % PLOT REFERENCE IF IT EXISTS
            try
                if isfield(options.mcmc,'m_ref');
                    hold on
                    try
                        plot(prior{im}.x,options.mcmc.m_ref{im},'r-','LineWidth',2);
                    catch                        
                        plot(options.mcmc.m_ref{im},'r-','LineWidth',2);
                        %    sippi_verbose(sprintf('cannot plot m_ref'));
                    end
                    hold off
                end
            catch
                keyboard
            end
            xlabel('X')
            ylabel(prior{im}.name)
            print_mul(sprintf('%s_m%d_posterior_sample_density',fname,im),options.plot.hardcopy_types)
            
            %%
            figure_focus(f_id);
            try
                plot(x,reals,'k-','linewidth',.1);
            catch
                plot(reals,'k-','linewidth',.1);
            end
            hold on
            %plot(prior{im}.x,etype_mean,'r-','linewidth',2);
            if exist('quantile','file')
                try
                    plot(x,quantile(reals',.025),'g--','linewidth',1);
                    plot(x,quantile(reals',.5),'g-','linewidth',2);
                    plot(x,quantile(reals',.975),'g--','linewidth',1);
                catch
                    plot(quantile(reals',.025),'g--','linewidth',1);
                    plot(quantile(reals',.5),'g-','linewidth',1);
                    plot(quantile(reals',.975),'g--','linewidth',1);
                end    
            end
            hold off
            xlabel('X')
            ylabel(prior{im}.name)
            % optonallt set y axis-limits
            if isfield(prior{im},'cax');
                set(gca,'ylim',prior{im}.cax);
            end
            
            % PLOT REFERENCE IF IT EXISTS
            try
                if isfield(options.mcmc,'m_ref');
                    hold on
                    try
                        plot(prior{im}.x,options.mcmc.m_ref{im},'r-','MarkerSize',11,'LineWidth',2);
                    catch
                        plot(options.mcmc.m_ref{im},'r-','MarkerSize',11,'LineWidth',2);
                        %sippi_verbose(sprintf('cannot plot m_ref'));
                    end
                    hold off
                end
            end
            
        else
            if ax_lscape==1;
                nsp_y=5;
                nsp_x=ceil(n_reals(im)/nsp_y);
            else
                nsp_x=5;
                nsp_y=ceil(n_reals(im)/nsp_x);
            end
            try;clear m;end
            for i=1:n_reals(im)
                
                subplot(nsp_y,nsp_x,i);
                
                use_colorbar=0;
                i_cb=ceil((nsp_y+1)/2)*nsp_x;
                if i==i_cb; use_colorbar=1;end
                
                try
                    if (length(z)>1)
                        m{i}{im}=reals(:,:,:,i);
                    else
                        m{i}{im}=reals(:,:,i);
                    end
                    sippi_plot_prior(prior,m{i},im,use_colorbar,f_id);
                catch
                    disp(sprintf('%s : failed to plot realization %d',mfilename,i))
                end
            end
        end
        if options.plot.suptitle==1,
            st=suptitle(title_txt);
            set(st,'interp','none','FontSize',18);
        else
            %title(title_txt)
        end
        print_mul(sprintf('%s_m%d',fname,im),options.plot.hardcopy_types)
        
        
        %% PLOT ETYPES
        if ndim>1
            f_id=(im-1)*10+2;
            figure_focus(f_id);set_paper('landscape');clf;
            
            % ETYPE MEAN
            if (ax_lscape==1)
                subplot(2,1,1);
            else
                subplot(1,2,1);
            end
            set(gca,'FontSize',options.plot.axis.fontsize)
            met{im}=etype_mean;
            sippi_plot_prior(prior,met,im,0,f_id);colorbar off;
            title('Posterior mean')
            try;caxis(cax);end
            cb=colorbar_shift;
            set(get(cb,'Ylabel'),'String','Sample Mean')
            
            % ETYPE VARIANCE
            if (ax_lscape==1)
                subplot(2,1,2);
            else
                subplot(1,2,2);
            end
            set(gca,'FontSize',options.plot.axis.fontsize)
            %met{im}=etype_var;
            met{im}=sqrt(etype_var);
            sippi_plot_prior(prior,met,im,0,f_id);colorbar off;
            title('Posterior std')
            xlabel('X');ylabel('Y');zlabel('Z')
            cax_var=[0 max(etype_var(:))];
            try
                Va=deformat_variogram(prior{im}.Va);
                cax_var(2)=sum(Va.par1);
                
            end
            %try;caxis(cax_var);end
            try;caxis(sqrt(cax_var));end
            cb=colorbar_shift;
            %set(get(cb,'Ylabel'),'String','Sample Variance')
            set(get(cb,'Ylabel'),'String','Sample std')
            % SUPTITLE
            if options.plot.suptitle==1,
                st=suptitle(sprintf('m%d: %s',im,prior{im}.name));
                set(st,'interpreter','none');
            end
            print_mul(sprintf('%s_m%d_sample_stat',fname,im),options.plot.hardcopy_types)
        end
        
        
        if ~exist('mcmc','var')
            cd(cwd)
            return
        end
        
        
        %% PLOT ACCEPTANCE RATE
        try
            
            fn=(im-1)*10+5;
            figure_focus(fn);set_paper('landscape');clf;
            
            % NEXT LINE ADDED TO ACCOUNT FOR MULTIPLE CHAINS
            mcmc=options.C{1}.mcmc;
            acc=mcmc.acc(im,1:mcmc.i);
            perturb=mcmc.perturb(im,1:mcmc.i);
            ip=find(perturb==1); % find indice of iteration when parameter has been perturbed
            
            fak=(1/10)*length(acc)./prior{im}.seq_gibbs.n_update_history;; % smoothing factor
            fak=1; % smoothing factor
            AccNum=conv_strip(acc(ip),ones(1,fak*prior{im}.seq_gibbs.n_update_history));
            AccRate_smooth=(1/fak)*AccNum/prior{im}.seq_gibbs.n_update_history;
            AccRate=AccNum/prior{im}.seq_gibbs.n_update_history;
            subplot(2,1,1);
            set(gca,'FontSize',options.plot.axis.fontsize)
            try;title(sprintf('m%d : %s',im,prior{im}.name));end
            plot(ip,AccRate_smooth,'-');
            xlabel('Iteration number');
            ylabel('Acceptance rate')
            title(title_txt)
            ylim=get(gca,'ylim');
            if ylim(2)>1.1; ylim(2)=1.1;end
            set(gca,'ylim',ylim);
            set(gca,'FontSize',options.plot.axis.fontsize)
            
            subplot(2,1,2);
            set(gca,'FontSize',options.plot.axis.fontsize)
            hist(AccRate,linspace(0,1,21));
            set(gca,'xlim',[0 1])
            xlabel('Acceptance Rate')
            ylabel('pdf')
            
            print_mul(sprintf('%s_m%d_rate',fname,im),options.plot.hardcopy_types)
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot acceptance rate for prior{%d}',mfilename,im));
            cd(cwd);
        end
        
        %% PLOT CORRELATION COEFFICIENT / FIND NITE PER INDEPENDANT POST REAL
        try
            
            if ndim==0
                %% autocorrelation analysis... to come
                fn=(im-1)*10+6;
                figure_focus(fn);set_paper('landscape');clf;
                set(gca,'FontSize',options.plot.axis.fontsize)
                
                c_reals=reals_all;
                
                % ONLY USE REALS AFTER SEQ GIBBS HAS FINISHED...
                try
                    c_i1=prior{im}.seq_gibbs.i_update_step_max/options.mcmc.i_sample;
                    c_reals=reals_all(c_i1:size(reals_all,1),:);
                end
                
                
                % compute cross correlation
                r=c_reals-mean(c_reals);
                c=conv(r,flip(r));
                
                c=c(length(c_reals):end);
                c=c./max(c);
                xc=[0:1:(length(c))-1].*options.mcmc.i_sample;
                plot(xc,c,'-','linewidth',2);grid on
                ic0=find(c<0);ic0=ic0(1);
                axis([0 xc(ic0)*8 -.5 1])
                hold on;
                plot([1 1].*xc(ic0),[-1 1]*.2,'-','linewidth',3);
                text(xc(ic0)+0.01*diff(get(gca,'xlim')),0.1,sprintf('Nite=%d',xc(ic0)),'FontSize',options.plot.axis.fontsize)
                hold off
                
                xlabel('iteration #')
                ylabel(sprintf('autocorrelation of %s(m%d)',prior{im}.name,im))
                set(gca,'FontSize',options.plot.axis.fontsize)
                
                print_mul(sprintf('%s_autocorr_m%d',fname,im),options.plot.hardcopy_types)
                
                %%
                
            elseif ndim>=1
                
                %% Ananlysis based on correlation coefficient
                fn=(im-1)*10+6;
                figure_focus(fn);set_paper('landscape');clf;
                set(gca,'FontSize',options.plot.axis.fontsize)
                nr=size(reals_all,1);
                it=[1:1:nr].*options.mcmc.i_sample;
                for i=1:nr;
                    c=corrcoef(reals_all(i,:),reals_all(nr,:));
                    cc(i)=c(2);
                end
                plot(it,cc,'k-','linewidth',2);
                xlabel('iteration')
                ylabel('Correlation coefficient')
                ylim=get(gca,'ylim');
                
                % FIND N_IT FOR INDEPENDANT REALS
                [hh,hx]=hist(cc,30);
                lev=hx(find(hh==max(hh)));
                i_threshold=max(find(cc<lev(1)));
                n_threshold=it(nr)-it(i_threshold);
                txt=sprintf('About %d iterations between independant realizations',n_threshold);
                t=text(.1,.9,txt,'units','normalized','FontSize',options.plot.axis.fontsize);
                i_independant=it(nr)-[n_threshold:n_threshold:it(nr)];
                try;set(gca,'xlim',[1 max(it)]);end
                for i=1:length(i_independant)
                    hold on
                    plot([1 1].*i_independant(i),ylim,'r-','LineWidth',1.5)
                    hold off
                end
                
                try;title(sprintf('m%d : %s',im,prior{im}.name),'interp','none');end
                print_mul(sprintf('%s_m%d_corrcoeff',fname,im),options.plot.hardcopy_types)
                
                
                %% ANALYSIS BASED ON ESS, TAU, AUTOCORRELATION
                try
                    fn=(im-1)*10+7;
                    figure_focus(fn);set_paper('landscape');clf;
                    %[e,ni]=ESS(reals_all(:,:)',400,1,options.mcmc.i_sample);
                    
                    n_max=1000;
                    di=max([ceil(size(reals_all,2)/n_max),1]);
                    [e,ni,e_av,ni_av]=ESS(reals_all(:,1:di:end),400,1,options.mcmc.i_sample*di);
                    disp(sprintf('%s: prior{%d}, Number of independent realizations: %d',mfilename,im,floor(min(e))))
                    disp(sprintf('%s: prior{%d}, tau (average): %5.1f',mfilename,im,e_av))
                    print_mul(sprintf('%s_m%d_autocorr_ess',fname,im),options.plot.hardcopy_types);
                    
                catch
                    disp(sptinf('%s: Could not perform ESS analysis',mfilename))
                end
            end
            
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot corrcoeff stats for prior{%d}',mfilename,im));
            cd(cwd);
        end % end plot correlection coefficient
        
    end
    cd(cwd);
    
%catch
%    keyboard
%    sippi_verbose(sprintf('%s : something went wrong',mfilename))
%    cd(cwd)
%end
cd(cwd)
