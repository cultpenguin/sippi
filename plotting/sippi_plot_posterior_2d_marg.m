% sippi_plot_posterior_2d_marg: plots 2D posterior marginal distributions
%
% Call:
%    [options,reals_all]=sippi_plot_posterior_2d_marg(options,prior,data,fname);
%
% See also: sippi_plot_posterior
%
function [options,reals_all]=sippi_plot_posterior_2d_marg(options,prior,data,fname);


%% LOAD THE CORRECT DATA
cwd=pwd;
if nargin==0
    % LOAD FROM MAT FILES
    [p,matfile]=fileparts(pwd);
    load(matfile);
    options.mcmc=mcmc;
elseif nargin==1;
    if isstruct(options),
    else
        fname=options;
        cd(fname);
        load(fname);
        options.mcmc=mcmc;
    end
else
    
    if nargin>3
        options.mcmc=mcmc;
    end
end
if nargin<5
    try
        fname=options.txt;
    catch
        fname=mfilename;
    end
end

% ALL DATA LOADED

%%
% SET DFAULT PLOTTING SETTINGS
options=sippi_plot_defaults(options);

pl_marg2d_scatter=options.plot.marg2d.pl_marg2d_scatter;
pl_marg2d_image=options.plot.marg2d.pl_marg2d_image;
pl_marg2d_scatter_combined=options.plot.marg2d.pl_marg2d_scatter_combined;
pl_marg2d_image_combined=options.plot.marg2d.pl_marg2d_image_combined;
pl_marg2d_hpd=options.plot.marg2d.pl_marg2d_hpd;

NX=options.plot.marg2d.NX;
NY=options.plot.marg2d.NY;

% find 1D type prior model
im_onedim=[];
for i=1:length(prior)
    if max(prior{i}.dim)==1
        im_onedim=[im_onedim, i];
    end
end

n=length(im_onedim);

if n==0
    pl_marg2d_scatter=0;
    pl_marg2d_image=0;
    pl_marg2d_scatter_combined=0;
    pl_marg2d_image_combined=0;
    pl_marg2d_hpd=0;
end

j=0;

%% GET REALS
clear reals*;
skip_seq_gibbs=options.plot.skip_seq_gibbs;
if isfield(options.mcmc,'adaptive_rejection');
    % use all realizations using sippi_rejection*
    skip_seq_gibbs=0;
end   
skip_seq_gibbs=1;

for k=1:(length(im_onedim))
    [reals1,etype_mean1,etype_var1,reals_all1]=sippi_get_sample('.',im_onedim(k),100000,skip_seq_gibbs);
    reals_all(:,k)=reals1(:);
    
    try
        cax(k,:)=[prior{k}.min prior{k}.max];
    catch
        cax(k,:)=[min(reals_all(:,k)) max(reals_all(:,k))];
    end
    
end
%%
for k=1:(length(im_onedim)-1)
    for l=(k+1):(length(im_onedim))
        
                
        if pl_marg2d_scatter==1;
            figure_focus(50+j);clf;set_paper('landscape');
            plot(reals_all(:,k),reals_all(:,l),'k.')
            scatter(reals_all(:,k),reals_all(:,l),1,1:size(reals_all,1))           
            try;xlabel(prior{im_onedim(k)}.name);end
            try;ylabel(prior{im_onedim(l)}.name);end
            
            % REF MODEL
            if isfield(options.mcmc,'m_ref');
                try
                    hold on;plot(options.mcmc.m_ref{k},options.mcmc.m_ref{l},'ro','MarkerSize',6,'LineWidth',3);hold off
                end
            end
            ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            set(gca,'xlim',cax(k,:));
            set(gca,'ylim',cax(l,:));
            print_mul(sprintf('%s_post_marg2d_m%d_m%d_scatter',fname,im_onedim(k),im_onedim(l)),options.plot.hardcopy_types)
            
        end
        
        %% 2d marg image
        if pl_marg2d_image==1;
            figure_focus(60+j);clf;set_paper('landscape');
            
            x_arr=linspace(cax(k,1),cax(k,2),NX);
            y_arr=linspace(cax(l,1),cax(l,2),NY);
            try
                % if prior{im}.min,prior{im}.max exists
                [Z,x_arr,y_arr] = hist2(reals_all(:,k),reals_all(:,l),x_arr,y_arr);
            catch
                [Z,x_arr,y_arr] = hist2(reals_all(:,k),reals_all(:,l),NX,NY);
            end
            
            if pl_marg2d_hpd==1
                levels=hpd_2d(Z,options.plot.marg2d.hpd_interval);
                Zl=Z.*0;
                for il=1:length(levels);
                    Zl(Z>levels(il))=il;
                end
                %contourf(Z,levels);
                imagesc(x_arr,y_arr,Zl');
            else
                imagesc(x_arr,y_arr,Z);
            end
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
            ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            set(gca,'xlim',cax(k,:));
            set(gca,'ylim',cax(l,:));
            
            if pl_marg2d_hpd==1;
                figure_name=sprintf('%s_post_marg2d_m%d_m%d_image_hpd',fname,im_onedim(k),im_onedim(l));
            else
                figure_name=sprintf('%s_post_marg2d_m%d_m%d_image',fname,im_onedim(k),im_onedim(l));
            end
            print_mul(sprintf(figure_name,fname),options.plot.hardcopy_types)
            
        end
    end
end

%% 2d marginals (image) on one plot
try
if (pl_marg2d_image_combined==1)&(length(im_onedim)>0);
    fn=figure_focus(71);clf;set_paper('landscape');
    for j=1:(n-1)
        for k=((1)+j):n
            
            r1=reals_all(:,j);
            r2=reals_all(:,k);
            x_arr=linspace(cax(j,1),cax(j,2),NX);
            y_arr=linspace(cax(k,1),cax(k,2),NY);
            
            try
                [Z,x_arr,y_arr] = hist2(r1(:),r2(:),x_arr,y_arr);
            catch                
                [Z,x_arr,y_arr] = hist2(r1(:),r2(:));
            end
            isp=(j-1)*(n-1)+(k-1);
            
            subplot(n-1,n-1,isp);
            if pl_marg2d_hpd==1
                levels=hpd_2d(Z,options.plot.marg2d.hpd_interval);
                Zl=Z.*0;
                for il=1:length(levels);
                    Zl(Z>levels(il))=il;
                end
                %contourf(Z,levels);
                imagesc(x_arr,y_arr,Zl');
            else
                imagesc(x_arr,y_arr,Z);
            end
            set(gca,'ydir','normal');
            %xlabel(prior{im_onedim(j)}.name,'interp','none')
            %ylabel(prior{im_onedim(k)}.name,'interp','none')
            xlabel(prior{im_onedim(j)}.name)
            ylabel(prior{im_onedim(k)}.name)
            
            set(gca,'xlim',cax(j,:));
            set(gca,'ylim',cax(k,:));
            
            
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
    if pl_marg2d_hpd==1;
        figure_name=sprintf('%s_post_marg2d_image_combined_hpd',fname);
    else
        figure_name=sprintf('%s_post_marg2d_image_combined',fname);
    end
    print_mul(sprintf(figure_name,fname),options.plot.hardcopy_types)
end
end
%% 2d marginals (scatter) on one plot
if (pl_marg2d_scatter_combined==1)&(length(im_onedim)>0);
    FS = max([12-n,6]);
    fn=figure_focus(72);clf;set_paper('landscape');
    for j=1:(n-1)
        for k=((1)+j):n
            
            r1=reals_all(:,j);
            r2=reals_all(:,k);
            isp=(j-1)*(n-1)+(k-1);
            subplot(n-1,n-1,isp);
            set(gca,'FontSize',FS)
            set(gca,'ydir','normal');
            %plot(reals_all(:,j),reals_all(:,k),'k.','MarkerSize',3)
            scatter(reals_all(:,j),reals_all(:,k),3,1:size(reals_all,1),'filled')           
            set(gca,'FontSize',FS)
            %xlabel(prior{im_onedim(j)}.name,'interp','none')
            %ylabel(prior{im_onedim(k)}.name,'interp','none')
            xlabel(prior{im_onedim(j)}.name,'FontSize',FS)
            ylabel(prior{im_onedim(k)}.name,'FontSize',FS)
            try
                if isfield(options.mcmc,'m_ref');
                    hold on
                    plot(options.mcmc.m_ref{j},options.mcmc.m_ref{k},'ro','MarkerSize',6,'LineWidth',3);
                    hold off
                end
            end
            box on
            try;set(gca,'xlim',cax(j,:));end
            try,set(gca,'ylim',cax(k,:));end
            %colormap(1-gray);
        end
    end
    print_mul(sprintf('%s_post_marg2d_scatter_combined',fname),options.plot.hardcopy_types)
end

%% GO BACK TO STARTING DIRECTORY
cd(cwd)
