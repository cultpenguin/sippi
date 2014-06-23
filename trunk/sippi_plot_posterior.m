function sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
% sippi_plot_posterior Plot statistics from posterior sample
%
% Call :
%    sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
%
% See also sippi_plot_prior
%

if nargin==0;
    [f1,fname]=fileparts(pwd);
end

if ~exist('supt','var');
    supt=0;
end



pl_logL=1;
pl_base=1;
pl_2d_marg=1;
pl_data=1;


% color codes
col=[
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    0 0 1
    .5 .5 .5
    ];


cwd=pwd;

%% DATA
if isstr(fname)
    try
        cd(fname);
        load([fname,'.mat']);
    catch
        load([fname,'.mat']);
    end
    
else
    data=fname;
    fname='lsq';
end

plotdir=pwd;
try
    fname=options.txt;
end

%if ~isfield(options,'FS')
%    options.FS=12;
%end

%% PERHAPS SET THE OPTIONS IN A SEPERATE MFILE
options.axis.null='';
if ~isfield(options.axis,'fontsize');options.axis.fontsize=20;end
if ~isfield(options.axis,'width');options.axis.width=8;end
if ~isfield(options.axis,'height');options.axis.height=8;end
if ~isfield(options.axis,'w0');options.axis.w0=2;end
if ~isfield(options.axis,'h0');options.axis.h0=2;end

prior=sippi_prior_init(prior);

%% logL curve
if pl_logL==1;
    
    figure(31);set_paper('landscape');
    set(gca,'FontSize',options.axis.fontsize);
    for i=1:length(prior);i_update_step_max(i)=prior{i}.seq_gibbs.i_update_step_max;end
    i_update_step_max=max(i_update_step_max);
    try
        sippi_plot_loglikelihood(mcmc.logL_all);
        legend(num2str([1:size(mcmc.logL_all,1)]'))
        y1=max(max(mcmc.logL_all));
        try
            y2=min(min(mcmc.logL_all(:,i_update_step_max)))
        catch
            y2=min(min(mcmc.logL_all));
        end
    catch
        sippi_plot_loglikelihood(mcmc.logL);
        y1=max(max(mcmc.logL));
        try
            y2=(min(mcmc.logL(:,i_update_step_max:end)));
        catch
            y2=(min(mcmc.logL));
        end
    end
    try
        xlim=get(gca,'xlim');
        % GET TOTAL NUMBER OF DATA
        N=0;for l=1:length(data);N=N+length(data{l}.d_obs);end;
        hold on
        plot(xlim,[1 1].*(-N/2),'r-')
        plot(xlim,[1 1].*(-N/2-2*sqrt(N/2)),'r--')
        plot(xlim,[1 1].*(-N/2+2*sqrt(N/2)),'r--')
        hold off
        if y2>(-N/2-4*sqrt(N/2)), y2=(-N/2-4.1*sqrt(N/2));end
        if y1<(-N/2+4*sqrt(N/2)), y1=(-N/2+4.1*sqrt(N/2));end
        
    end
    try
        set(gca,'ylim',[y2 y1])
    end
     set(gca,'FontSize',options.axis.fontsize)           
    print_mul(sprintf('%s_logL',fname))
    xlim=get(gca,'xlim');
    set(gca,'xlim',[xlim(1) xlim(2)/20]);
    print_mul(sprintf('%s_logL_start',fname))
    
    
    %% autocorrelation
    try
        figure(3);clf;set_paper('landscape');
        
        %nite=length(mcmc.logL);
        i1=1;
        for i=1:length(prior);
            i1 = max([prior{i}.seq_gibbs.i_update_step_max i1]);
        end
        
        ii=i1:length(mcmc.logL);
        c=xcorr(mcmc.logL(ii)-mean(mcmc.logL(ii)));
        c=c(length(ii):end);
        c=c./max(c);
        xc=[0:1:(length(c))-1];
        plot(xc,c,'-');grid on
        ic0=find(c<0);ic0=ic0(1);
        axis([0 xc(ic0)*8 -.5 1])
        hold on;
        plot([1 1].*xc(ic0),[-1 1]*.2,'-','linewidth',3);
        text(xc(ic0)+0.01*diff(get(gca,'xlim')),0.1,sprintf('Nite=%d',xc(ic0)),'FontSize',options.axis.fontsize)
        hold off
        
        xlabel('iteration #')
        ylabel('autocorrelation of logL')
        set(gca,'FontSize',options.axis.fontsize)           
        print_mul(sprintf('%s_logL_autocorr',fname))
        
    end
end



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
            n_reals(j)=15;
        end
    end
end

if length(n_reals)==1;
    n_reals=ones(1,length(prior)).*n_reals;
end



%for im=im_arr;
if pl_base==1;
    
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
        
        
        [reals,etype_mean,etype_var,reals_all,ite_reals]=sippi_get_sample(data,prior,id,im,n_reals(im),options);
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
                ii=[1:1:length(reals_all)].*options.mcmc.i_sample;
                plot(ii(:),reals_all);
                
                ylim=get(gca,'ylim');
                if isfield(prior{im},'cax');ylim=prior{im}.cax;end
                if isfield(prior{im},'min');ylim(1)=prior{im}.min;end
                if isfield(prior{im},'max');ylim(2)=prior{im}.max;end
                set(gca,'ylim',ylim);
                set(gca,'FontSize',options.axis.fontsize)
                xlabel('Iteration #')
                ylabel(prior{im}.name)
                print_mul(sprintf('%s_m%d_posterior_values',fname,im))
                
            catch
                disp(sprintf('%s : failed to plot 1D posterior values for prior #%d',mfilename,im));
            end
            
            figure_focus(f_id);
            %% continue
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
            end
            
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
            xlabel(prior{im}.name,'interpreter','none','FontSize',options.axis.fontsize+2)
            ylabel('Frequency','interpreter','none','FontSize',options.axis.fontsize+2)
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
            
            %ppp(options.axis.width,options.axis.height,options.axis.fontsize,options.axis.w0,options.axis.h0);
            
        elseif ndim==1
            plot(prior{im}.x,reals,'k-');
            hold on
            %plot(prior{im}.x,etype_mean,'r-','linewidth',2);
            plot(prior{im}.x,quantile(reals',.025),'r--','linewidth',2);
            plot(prior{im}.x,quantile(reals',.5),'r-','linewidth',2);
            plot(prior{im}.x,quantile(reals',.975),'r--','linewidth',2);
            hold off
            xlabel('X')
            ylabel(prior{im}.name)
            % optonallt set y axis-limits
            if isfield(prior{im},'cax');
                set(gca,'ylim',prior{im}.cax);
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
                
                progress_txt(i,n_reals(im),'computing data response')
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
                    %sippi_plot_prior(prior,m{i}{id},im,use_colorbar,f_id);
                catch
                    disp(sprintf('%s : failed to plot realization %d',mfilename,i))
                end
            end
        end
        if supt==1,
            st=suptitle(title_txt);
            set(st,'interp','none','FontSize',18);
        else
            %title(title_txt)
        end
        print_mul(sprintf('%s_m%d_posterior_sample',fname,im))
        
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
            set(gca,'FontSize',options.axis.fontsize)
            met{im}=etype_mean;
            sippi_plot_prior(prior,met,im,0,f_id);colorbar off;
            caxis(cax);
            cb=colorbar_shift;
            set(get(cb,'Ylabel'),'String','Sample Mean')
            
            % ETYPE VARIANCE
            if (ax_lscape==1)
                subplot(2,1,2);
            else
                subplot(1,2,2);
            end
            set(gca,'FontSize',options.axis.fontsize)
            %met{im}=etype_var;
            met{im}=sqrt(etype_var);
            sippi_plot_prior(prior,met,im,0,f_id);colorbar off;
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
            set(get(cb,'Ylabel'),'String','Sample STD')
            
            % SUPTITLE
            if supt==1,
                st=suptitle(sprintf('m%d: %s',im,prior{im}.name));
                set(st,'interpreter','none');
            end
            print_mul(sprintf('%s_m%d_sample_stat',fname,im))
        end
        
        
        if ~exist('mcmc','var')
            cd(cwd)
            
            return
        end
        
        
        %% PLOT ACCEPTANCE RATE
        try
            
            fn=(im-1)*10+5;
            figure_focus(fn);set_paper('landscape');clf;
            acc=mcmc.acc(im,1:mcmc.i);
            perturb=mcmc.perturb(im,1:mcmc.i);
            ip=find(perturb==1); % find indice of iteration when parameter has been perturbed
            
            fak=(1/10)*length(acc)./prior{im}.seq_gibbs.n_update_history;; % smoothing factor
            fak=1; % smoothing factor
            AccNum=conv_strip(acc(ip),ones(1,fak*prior{im}.seq_gibbs.n_update_history));
            AccRate_smooth=(1/fak)*AccNum/prior{im}.seq_gibbs.n_update_history;
            AccRate=AccNum/prior{im}.seq_gibbs.n_update_history;
            subplot(2,1,1);
            try;title(sprintf('m%d : %s',im,prior{im}.name));end
            plot(ip,AccRate_smooth,'-');
            xlabel('Iteration number');
            ylabel('Acceptance rate')
            title(title_txt)
            ylim=get(gca,'ylim');
            if ylim(2)>1.1; ylim(2)=1.1;end
            set(gca,'ylim',ylim);
            subplot(2,1,2);
            hist(AccRate,linspace(0,1,21));
            set(gca,'xlim',[0 1])
            xlabel('Acceptance Rate')
            ylabel('pdf')
            
            print_mul(sprintf('%s_m%d__rate',fname,im))
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot acceptance rate',mfilename));
            cd(cwd);
        end
        %% PLOT CORRELATION COEFFICIENT / FIND NITE PER INDEPENDANT POST REAL
        try
            if ndim==0
                %% autocorrelation analysis... to come
                fn=(im-1)*10+6;
                figure_focus(fn);set_paper('landscape');clf;
                set(gca,'FontSize',options.axis.fontsize)
                
                c_reals=reals_all;
                
                % ONLY USE REALS AFTER SEQ GIBBS HAS FINISHED...
                try
                    c_i1=prior{im}.seq_gibbs.i_update_step_max/options.mcmc.i_sample;
                    c_reals=reals_all(c_i1:size(reals_all,1),:);
                end
                
                c=xcorr(c_reals-mean(c_reals));
                c=c(length(c_reals):end);
                c=c./max(c);
                xc=[0:1:(length(c))-1].*options.mcmc.i_sample;
                plot(xc,c,'-','linewidth',2);grid on
                ic0=find(c<0);ic0=ic0(1);
                axis([0 xc(ic0)*8 -.5 1])
                hold on;
                plot([1 1].*xc(ic0),[-1 1]*.2,'-','linewidth',3);
                text(xc(ic0)+0.01*diff(get(gca,'xlim')),0.1,sprintf('Nite=%d',xc(ic0)),'FontSize',options.axis.fontsize)
                hold off
                
                xlabel('iteration #')
                ylabel(sprintf('autocorrelation of %s(m%d)',prior{im}.name,im))
                set(gca,'FontSize',options.axis.fontsize)           
        
                print_mul(sprintf('%s_autocorr_m%d',fname,im))
                
                %%
                
            elseif ndim>1
                fn=(im-1)*10+6;
                figure_focus(fn);set_paper('landscape');clf;
                set(gca,'FontSize',options.axis.fontsize)
                nr=size(reals_all,1);
                it=[1:1:nr].*mcmc.i_sample;
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
                t=text(.1,.9,txt,'units','normalized','FontSize',options.axis.fontsize);
                i_independant=it(nr)-[n_threshold:n_threshold:it(nr)];
                try;set(gca,'xlim',[1 options.mcmc.nite]);end
                for i=1:length(i_independant)
                    hold on
                    plot([1 1].*i_independant(i),ylim,'r-','LineWidth',1.5)
                    hold off
                end
                
                try;title(sprintf('m%d : %s',im,prior{im}.name),'interp','none');end
                print_mul(sprintf('%s_m%d_corrcoeff',fname,im))
                
            end
            
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot corrcoeff stats',mfilename));
            cd(cwd);
        end
        
        
        %% PLOT  LOGL
        try
            if nm==1
                fn=(im-1)*10+9;
                figure_focus(fn);set_paper('landscape');clf;
                set(gca,'FontSize',options.axis.fontsize);
                sippi_plot_loglikelihood(mcmc.logL(1:mcmc.i),mcmc.acc(im,1:mcmc.i));
                smcmc=sort(mcmc.logL(1:mcmc.i));y_min=smcmc(ceil(mcmc.i/200));
                ylim=get(gca,'ylim');
                ylim(1)=y_min;
                ylim(2)=max(mcmc.logL).*.8;
                set(gca,'ylim',ylim),
                grid on
                title(options.txt,'interpreter','none');
                print_mul(sprintf('%s_logL',fname))
                %else
                %    figure_focus((im-1)*10+9);set_paper('landscape');clf;
                %    title(sprintf('m%d : %s',im,prior{im}.name))
            end
        catch
            try;close(fn);end
            disp(sprintf('%s : could not plot logL curve',mfilename));
            cd(cwd);
        end
        
    end
end

%% 2D POSTERIOR MARGINALS.
if (length(prior)<2); pl_2d_marg=0;end
if (pl_2d_marg==1),
    %try
        try;cd(plotdir);end
        im_onedim=[];
        for j=1:length(im_arr);
            if max(prior{j}.dim)==1
                im_onedim=[im_onedim, j];
            end
        end
        n=length(im_onedim);
        j=0;
        clear reals*;
        for k=1:(length(im_onedim)-1)
            [reals1,etype_mean1,etype_var1,reals_all1]=sippi_get_sample(data,prior,id,im_onedim(k),100000,options);
            %reals_all(:,k)=reals_all1(:);
            reals_all(:,k)=reals1(:);
            for l=(k+1):(length(im_onedim))
                j=j+1;
                
                %% 2d marg scatter
                [reals2,etype_mean2,etype_var2,reals_all2]=sippi_get_sample(data,prior,id,im_onedim(l),100000,options);
                %reals_all(:,l)=reals_all2(:);
                reals_all(:,l)=reals2(:);
                
                figure_focus(50+j);clf;set_paper('landscape');
                plot(reals_all(:,k),reals_all(:,l),'k.')
                try;xlabel(prior{im_onedim(k)}.name);end
                try;ylabel(prior{im_onedim(l)}.name);end
                
                try;set(gca,'xlim',[prior{im_onedim(k)}.min prior{im_onedim(k)}.max]);end
                try;set(gca,'ylim',[prior{im_onedim(l)}.min prior{im_onedim(l)}.max]);end
                % REF MODEL
                if isfield(options.mcmc,'m_ref');
                    try
                        hold on;plot(options.mcmc.m_ref{k},options.mcmc.m_ref{l},'ro','MarkerSize',6,'LineWidth',3);hold off
                    end
                end
                ppp(options.axis.width,options.axis.height,options.axis.fontsize,options.axis.w0,options.axis.h0);
                print_mul(sprintf('%s_post_marg_m%d_m%d_scatter',fname,im_onedim(k),im_onedim(l)));
                %% 2d marg image
                pl_2d_marg_image=0;
                if pl_2d_marg_image==1;
                    figure_focus(60+j);clf;set_paper('landscape');
                    try;
                        NX=ceil(sqrt(length(reals1)));
                        %NX=40;
                        NX=21;
                        NY=NX;
                        try
                            % if prior{im}.min,prior{im}.max exists
                            [Z,x_arr,y_arr] = hist2(reals_all(:,k),reals_all(:,l),linspace(prior{im_onedim(k)}.min,prior{im_onedim(k)}.max,NX),linspace(prior{im_onedim(l)}.min,prior{im_onedim(l)}.max,NY));
                        catch
                            [Z,x_arr,y_arr] = hist2(reals_all(:,k),reals_all(:,l),NX,NY);
                        end
                    catch
                        [Z,x_arr,y_arr] = hist2(reals1',reals2');
                    end
                    
                    imagesc(x_arr,y_arr,Z');
                    try;xlabel(prior{im_onedim(k)}.name);end
                    try ylabel(prior{im_onedim(l)}.name);end
                    % REF MODEL
                    if isfield(options.mcmc,'m_ref');
                        try
                            hold on;plot(options.mcmc.m_ref{k},options.mcmc.m_ref{l},'ro','MarkerSize',6,'LineWidth',3);hold off
                        end
                    end
                    
                    colormap(1-gray);
                    set(gca,'ydir','normal');
                    %colorbar
                    ppp(options.axis.width,options.axis.height,options.axis.fontsize,options.axis.w0,options.axis.h0);
                    print_mul(sprintf('%s_post_marg_hist_m%d_m%d_gray',fname,im_onedim(k),im_onedim(l)))
                end
            end
        end
        
        %% 2d marginals on one plot
        fn=figure_focus(70);clf;set_paper('landscape');
        for j=1:(n-1)
            for k=((1)+j):n
                
                r1=reals_all(:,j);
                r2=reals_all(:,k);
                try
                    NX=25;
                    NY=NX;
                    try
                        % if prior{im}.min,prior{im}.max exists
                        [Z,x_arr,y_arr] = hist2(r1(:),r2(:),linspace(prior{im_onedim(j)}.min,prior{im_onedim(j)}.max,NX),linspace(prior{im_onedim(k)}.min,prior{im_onedim(k)}.max,NY));
                    catch
                        [Z,x_arr,y_arr] = hist2(r1(:),r2(:),NX,NY);
                    end
                catch
                    [Z,x_arr,y_arr] = hist2(r1(:),r2(:));
                end
                levels=hpd_2d(Z,[.1,.5,.9]);
                Zl=Z.*0;
                for il=1:length(levels);
                    Zl(Z>levels(il))=il;
                end
                %contourf(Z,levels);
                isp=(j-1)*(n-1)+(k-1);
                subplot(n-1,n-1,isp);
                imagesc(x_arr,y_arr,Zl');
                set(gca,'ydir','normal');
                %plot(reals_all(:,j),reals_all(:,k),'k.','MarkerSize',.01)
                xlabel(prior{im_onedim(j)}.name,'interp','none')
                ylabel(prior{im_onedim(k)}.name,'interp','none')
                
                try
                    if isfield(options.mcmc,'m_ref');
                        hold on
                        plot(options.mcmc.m_ref{j},options.mcmc.m_ref{k},'ro','MarkerSize',6,'LineWidth',3);
                        hold off
                    end
                end
                
                colormap(1-gray);
            end
        end
        print_mul(sprintf('%s_post_marg_hist',fname))
    %catch
%         try;close(fn);end
%         fprintf('%s : could not plot 2D marginals\n',mfilename);
%         cd(cwd);
%         keyboard
    %end
end

%% PLOT DATA ASSOCIATED TO REALS
if pl_data==1,
    try
        
        %% THIS ONE NEADS SOME HEAVY EDITING TO HANDLE TWO DATA SETS!!!
        try;cd(plotdir);end
        
        %%
        nd=length(data);
        
        %% MFILE TO OBTAIN m_sample_1, m_real_2
        clear m_post
        N=15;
        for im=1:length(prior);
            [reals]=sippi_get_sample(data,prior,id,im,N+1,options);
            for j=1:N;
                m_post{j}{im}=reals(:,j+1);
            end
        end
        for j=1:N
            [d_real{j}]=sippi_forward(m_post{j},forward,prior,data);
        end
        
        %%
        
        for id=1:nd;
            
            %% PLOT REALIZATOIN OF NOISE
            if isfield(data{id},'d0');
                noise_real=gaussian_simulation_cholesky(data{1}.d0,data{1}.CD,N);
            else
                noise_real=gaussian_simulation_cholesky(0,data{1}.CD,N);
            end
            
            %% GET DATA AND RESIDUAL FOR 5 REALIZATIONS FROM POST
            for i=1:N;
                data_obs(:,i)=data{1}.d_obs;
                data_forward(:,i)=d_real{i}{id};
                data_res(:,i)=data{1}.d_obs-d_real{i}{id};
            end
            
            %%
            f_handle=70+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            wiggle(1:N,1:size(noise_real,1),data_forward,'wiggle',.1);
            hold on;
            wiggle(1:N,1:size(noise_real,1),data_obs,'wiggle',.1);
            hold off
            set(gca,'FontSize',options.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            print_mul(sprintf('%s_id%d_post',fname,id))
            
            %%
            f_handle=80+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),data_res,'VA',.1);
            set(gca,'FontSize',options.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Data residual')
            print_mul(sprintf('%s_id%d_res',fname,id))
            %%
            f_handle=85+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),noise_real,'VA',.1);
            set(gca,'FontSize',options.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Noise realization')
            print_mul(sprintf('%s_id%d_noise',fname,id))
            %%
            f_handle=90+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),noise_real,'wiggle',.1);
            hold on
            wiggle(1:N,1:size(noise_real,1),data_res,'wiggle',.1);
            hold off
            set(gca,'FontSize',options.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Realizations of noise model (black) and posterior data residuals (red)')
            print_mul(sprintf('%s_id%d_res_noise',fname,id))
            
            %%
            f_handle=95+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            plot(data{id}.CD(:,1),'k-','LineWidth',2)
            set(gca,'xlim',[0 size(data{id}.CD,1)/4])
            xlabel('data #')
            ylabel('Covariance')
            %title('')
            ppp(options.axis.width,options.axis.height,options.axis.fontsize,options.axis.w0,options.axis.h0);
            print_mul(sprintf('%s_id%d_CD',fname,id))
            
            %%
            f_handle=95+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            bar(1:N,var(data_res));
            hold on
            plot([1 N],[1 1].*data{1}.CD(1),'k:','LineWidth',12)
            hold off
            set(gca,'xlim',[0 N+1])
            xlabel('realization #')
            ylabel('variance')
            %title('')
            ppp(options.axis.width,options.axis.height,options.axis.fontsize,options.axis.w0,options.axis.h0);
            print_mul(sprintf('%s_id%d_var_check',fname,id))
            
            
            %%
        end
        
    catch
        keyboard
        cd(cwd);
        close(f_handle)
        fprintf('%s : Cannot plot data response. \n',mfilename)
    end
end



%%
cd(cwd);



