function sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
% sippi_plot_posterior Plot statistics from posterior sample
%
% Call :
%    sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
%
% See also sippi_plot_prior
%



if ~exist('supt','var');
    supt=0;
end


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

if ~isfield(options,'FS')
    options.FS=12;
end

%%
options.axis_fontsize=8;
options.width=10;
options.height=10;
options.w0=2;
options.h0=2;


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



for im=im_arr;
    
    try;cd(plotdir);end
    clear cax;
    % find dimension
    ndim=sum(prior{im}.dim>1);
    %if ndim==0;ndim=1;end
    
    
    options.null='';
    id=1;
    
    [reals,etype_mean,etype_var,reals_all]=sippi_get_sample(data,prior,id,im,n_reals(im),options);
    m_post{im}=reals;
    
    if ~exist('cax','var');
        if isfield(prior{im},'cax')
            cax=prior{im}.cax;
        else
            try
                cax=[prior{im}.min prior{im}.max];
            catch
                cax=[min(reals(:)) max(reals(:))];
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
        N=length(reals);
        %N=n_reals(im);
        prior_sample=zeros(1,N);
        clear p;
        p{1}=prior{im};
        for i=1:n_reals(im);
            m=sippi_prior(p);
            sample_prior(i)=m{1};
        end
        
        hx=linspace(cax(1),cax(2),30);
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
        xlabel(prior{im}.name,'interpreter','none')
        legend('prior','posterior')
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
     
        
        ppp(options.width,options.height,options.axis_fontsize,options.w0,options.h0);
        
        %set(gca,'FontSize',16),
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
        if (prior{im}.dim(1)>max(prior{im}.dim(2:3)))
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
            if ((n_reals==i)&(i==(1*nsp_x)))|(i==(2*nsp_x));
                %if (i==(2*nsp_x));
                use_colorbar=1;
            end
            
            try
                if (length(z)>1)
                    m{i}{im}=reals(:,:,:,i);
                else
                    m{i}{im}=reals(:,:,i);
                end
                sippi_plot_model(prior,m{i},im,use_colorbar,f_id);
                %sippi_plot_model(prior,m{i}{id},im,use_colorbar,f_id);
            catch
                disp(sprintf('%s : failed to plot realization %d',mfilename,i))
            end
        end
    end
    if supt==1,
        st=suptitle(sprintf('m%d : %s',im,prior{im}.name));
        set(st,'interp','none','FontSize',18);
    else
        %title(sprintf('m#%d, posterior',im))
    end
    print_mul(sprintf('%s_m%d_posterior_sample',fname,im))
    
    
    %% PLOT ETYPES
    if ndim>1
        f_id=(im-1)*10+2;
        figure_focus(f_id);set_paper('landscape');clf;
        if (prior{im}.dim(1)>max(prior{im}.dim(2:3)))
            subplot(2,1,1);
        else
            subplot(1,2,1);
        end
        set(gca,'FontSize',options.FS)
        met{im}=etype_mean;
        sippi_plot_model(prior,met,im,0,f_id);colorbar off;
        caxis(cax);
        colorbar_shift;
        axis image
        title('sample mean')
        
        if (prior{im}.dim(1)>max(prior{im}.dim(2:3)))
            subplot(2,1,2);
        else
            subplot(1,2,2);
        end
        set(gca,'FontSize',options.FS)
        met{im}=etype_var;
        sippi_plot_model(prior,met,im,0,f_id);colorbar off;
        xlabel('X');ylabel('Y');zlabel('Z')
        cax_var=[0 max(etype_var(:))];
        try
            Va=deformat_variogram(prior{im}.Va);
            cax_var(2)=sum(Va.par1);
            
        end
        try;caxis(cax_var);end
        colorbar_shift;
        axis image
        
        title('sample variance')
        if supt==1,
            st=suptitle(sprintf('m%d : %s',im,prior{im}.name));
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
        title(sprintf('m#%d, posterior',im))
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
        
        if ndim>1
            fn=(im-1)*10+6;
            figure_focus(fn);set_paper('landscape');clf;
            set(gca,'FontSize',options.FS)
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
            t=text(.1,.9,txt,'units','normalized','FontSize',options.FS);
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
            set(gca,'FontSize',options.FS);
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


%% 2D POSTERIOR MARGINALS.
try
    im_onedim=[];
    for j=1:length(im_arr);
        if max(prior{j}.dim)==1
            im_onedim=[im_onedim, j];
        end
    end
    n=length(im_onedim);
    for k=1:(length(im_onedim)-1)
        [reals1,etype_mean1,etype_var1,reals_all1]=sippi_get_sample(data,prior,id,im_onedim(k),100000,options);
        [reals2,etype_mean2,etype_var2,reals_all2]=sippi_get_sample(data,prior,id,im_onedim(k+1),100000,options);
        
        reals_all(:,k)=reals_all1(:);
        reals_all(:,k+1)=reals_all2(:);
        
        %keyboard
        figure_focus(50+k);clf;set_paper('landscape');
        plot(reals_all1,reals_all2,'.')
        try;xlabel(prior{im_onedim(k)}.name);end
        try;ylabel(prior{im_onedim(k+1)}.name);end
        
        try;set(gca,'xlim',[prior{im_onedim(k)}.min prior{im_onedim(k)}.max]);end
        try;set(gca,'ylim',[prior{im_onedim(k+1)}.min prior{im_onedim(k+1)}.max]);end
        ppp(options.width,options.height,options.axis_fontsize,options.w0,options.h0);
        print_mul(sprintf('%s_post_marg_m%d_m%d',fname,im_onedim(k),im_onedim(k+1)));
        
        figure_focus(60+k);clf;set_paper('landscape');
        try;
            NX=ceil(sqrt(length(reals1)));
            %NX=40;
            NX=21;
            NY=NX;
            try
                % if prior{im}.min,prior{im}.max exists
                [Z,x_arr,y_arr] = hist2(reals_all1(:),reals_all2(:),linspace(prior{im_onedim(k)}.min,prior{im_onedim(k)}.max,NX),linspace(prior{im_onedim(k+1)}.min,prior{im_onedim(k+1)}.max,NY));
            catch
                [Z,x_arr,y_arr] = hist2(reals_all1(:),reals_all2(:),NX,NY);
            end
        catch
            [Z,x_arr,y_arr] = hist2(reals1',reals2');
        end
        
        imagesc(x_arr,y_arr,Z');
        try;xlabel(prior{im_onedim(k)}.name);end
        try ylabel(prior{im_onedim(k+1)}.name);end
        
        try
            if isfield(options.mcmc,'m_ref');
                hold on
                plot(options.mcmc.m_ref{k},options.mcmc.m_ref{k+1},'ro','MarkerSize',6,'LineWidth',3);
                hold off
            end
        end
          
        
        colormap(1-gray);
        set(gca,'ydir','normal');
        %colorbar
        ppp(options.width,options.height,options.axis_fontsize,options.w0,options.h0);
        print_mul(sprintf('%s_post_marg_hist_m%d_m%d',fname,im_onedim(k),im_onedim(k+1)))
        
    end
    
    
    %% 2d marginals on one plot
    figure_focus(70);clf;set_paper('landscape');
    for j=1:(n-1)
        for k=((1)+j):n
            
            r1=reals_all(:,j);r2=reals_all(:,k);
            try
                %NX=ceil(1*sqrt(length(reals1)));
                %NX=40;
                NX=15;
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
catch
    try;close(fn);end
    fprintf('%s : could not plot 2D marginals\n',mfilename);
    cd(cwd);
    keyboard
end

%% PLOT DATA ASSOCIATED TO REALS
try
    
    %% THIS ONE NEADS SOME HEAVY EDITING TO HANDLE TWO DATA SETS!!!
    try;cd(plotdir);end
    
    %%
    f_handle=(im-1)*10+3;
    figure_focus(f_handle);set_paper('landscape');clf;
    subplot(1,1,1);
    set(gca,'FontSize',options.FS);
    nd=length(data);
    for id=1:nd;
        %if ~isfield(data{id},'i_use'); data{id}.i_use=1:1:(prod(size(data{id}.d_obs)));end
        subplot(1,nd,id)
        
        np=size(m_post{im},2);
        ii=ceil(linspace(1,np,min([50 np])));
        for k=ii:n_reals;
            [d_real,forward]=sippi_forward(m{k},forward,prior,data,id);
            p(1)=plot(d_real{id}(:),'-','col',col(1,:));
            d_obs=data{id}.d_obs(data{id}.i_use);
            dd(k,:)=d_obs(:)-d_real{id}(:);
            hold on
        end
        p(2)=plot(data{id}.d_obs(data{id}.i_use),'-*','col',col(2,:),'MarkerSize',2);
        hold off
        set(gca,'ylim',[min(data{id}.d_obs(:)).*.95 max(data{id}.d_obs(:)).*1.05])
        title(sprintf('Data #%d',id))
    end
    legend([p(1) p(2)],'d_{post}','d_{obs}')
    print_mul(sprintf('%s_d%d',fname,id))
    
    f_handle=(im-1)*10+4;
    figure_focus(f_handle);set_paper('landscape');clf;
    set(gca,'FontSize',options.FS);
    hist(dd(:),30);
    colormap(gray)
    xlabel('d_{obs}-d_{post}')
    ylabel('pdf')
    print_mul(sprintf('%s_d%d_posterior_datafit_hist',fname,id))
    %%
    
catch
    cd(cwd);
    close(f_handle)
    fprintf('%s : Cannot plot data response. \n',mfilename)
end




%%
cd(cwd);

return

%%
doSimCD=1;
if doSimCD==1;
    if strcmp(data.noise_model,'gaussian')
        clear L
        nsim=5000;
        cd_sim=gaussian_simulation_cholesky(0,data.CD,nsim);
        for k=1:nsim
            dd=cd_sim(data.i_use,k);
            logL(k)=-.5*dd'*data.iCD*dd;
        end
        hold on
        plot([xlim],[1 1].*mean(logL),'r-','LineWidth',2);
        plot([xlim],[1 1].*mean(logL)-2*std(logL),'r-','LineWidth',1);
        plot([xlim],[1 1].*mean(logL)+2*std(logL),'r-','LineWidth',1);
        hold off
        
        ylim=get(gca,'ylim');
        ylim(1)=max([ylim(1) mean(logL)-2*std(logL)]);
        try
            ylim(1)=min([ylim(1) 1.5*logL_Ref]);
        end
        set(gca,'ylim',ylim);
        
    end
end
ppp(13,10,12);
print_mul(sprintf('%s_logl_curve',fname))

cd(cwd);


