% sippi_plot_posterior_data: plots posterior data and noise relaizations
%
% Call 
%    [options]=sippi_plot_posterior_data(options,prior,data,forward);
%       
% See also: sippi_plot_posterior
%
function [options]=sippi_plot_posterior_data(options,prior,data,forward);


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
else
    
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

pl_obs_forward=1;
pl_data_res=1;
pl_noise_real=1;
pl_noise_data_real=1;
pl_cd=1;
pl_cd2=1;
pl_hist=1;


%% THIS ONE NEADS SOME HEAVY EDITING TO HANDLE TWO DATA SETS!!!
try;cd(plotdir);end

%%
nd=length(data);

try
    
    for id=1:nd;
        
        clear m_post
        N=15;
        for im=1:length(prior);
            [reals]=sippi_get_sample(data,prior,id,im,N+1,options);
            for j=1:N;
                if ndims(reals)==2;
                    m_post{j}{im}=reals(:,j+1);
                elseif ndims(reals)==3;
                    m_post{j}{im}=reals(:,:,j+1);
                else
                    m_post{j}{im}=reals(:,:,:,j+1);
                end
            end
        end
        for j=1:N
            [d_real{j}]=sippi_forward(m_post{j},forward,prior,data);
        end
        
        
        %% PLOT REALIZATOIN OF NOISE
        if isfield(data{id},'d0');
            noise_real=gaussian_simulation_cholesky(data{1}.d0,data{1}.CD,N);
        else
            noise_real=gaussian_simulation_cholesky(0,data{1}.CD,N);
        end
        noise_real=noise_real(data{id}.i_use,:);
        
        %% GET DATA AND RESIDUAL FOR 5 REALIZATIONS FROM POST
        for i=1:N;
            data_obs(:,i)=data{id}.d_obs(data{1}.i_use);
            data_forward(:,i)=d_real{i}{id};
            data_res(:,i)=data{id}.d_obs(data{1}.i_use)-d_real{i}{id};
        end
        
        %%
        if pl_obs_forward==1;
            f_handle=70+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            %wiggle(1:N,1:size(noise_real,1),data_forward,'wiggle',.1);
            plot(data_forward,'k-');
            hold on;
            %wiggle(1:N,1:size(noise_real,1),data_obs,'wiggle',.1);
            plot(data_obs,'r-');
            set(gca,'xlim',[1 size(data_obs,1)]);
            hold off
            set(gca,'FontSize',options.plot.axis.fontsize);
            xlabel('#')
            ylabel(sprintf('data #%d',id))
            print_mul(sprintf('%s_id%d_post',fname,id))
        end
        if pl_data_res==1;
            f_handle=80+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),data_res);
            set(gca,'FontSize',options.plot.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Data residual')
            print_mul(sprintf('%s_id%d_res',fname,id))
        end
        %%
        if pl_noise_real==1;
            f_handle=85+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),noise_real);
            set(gca,'FontSize',options.plot.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Noise realization')
            print_mul(sprintf('%s_id%d_noise',fname,id))
        end
        %%
        if pl_noise_data_real==1;
            f_handle=90+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            subplot(1,1,1);
            wiggle(1:N,1:size(noise_real,1),noise_real);
            hold on
            wiggle(1:N,1:size(noise_real,1),data_res);
            hold off
            set(gca,'FontSize',options.plot.axis.fontsize);
            xlabel('realization #')
            ylabel('data #')
            title('Realizations of noise model (black) and posterior data residuals (red)')
            print_mul(sprintf('%s_id%d_res_noise',fname,id))
        end
        %%
        if pl_cd==1;
            f_handle=95+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            plot(data{id}.CD(:,1),'k-','LineWidth',2)
            set(gca,'xlim',[0 size(data{id}.CD,1)/4])
            xlabel('data #')
            ylabel('Covariance')
            %title('')
            ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            print_mul(sprintf('%s_id%d_CD',fname,id))
        end
        %%
        if pl_cd2==1;
            f_handle=96+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            bar(1:N,var(data_res));
            hold on
            plot([1 N],[1 1].*data{1}.CD(1),'k:','LineWidth',12)
            hold off
            set(gca,'xlim',[0 N+1])
            xlabel('realization #')
            ylabel('variance')
            %title('')
            ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            print_mul(sprintf('%s_id%d_var_check',fname,id))
        end
        if pl_hist==1;
            f_handle=97+id;
            figure_focus(f_handle);set_paper('landscape');clf;
            [h_c,x_n]=hist(noise_real(:),30);
            [res_c]=hist(data_res(:),30);
            plot(x_n,h_c);
            hold on;
            plot(x_n,res_c,'k-','LineWidth',2)
            hold off
            legend('Noise realization','Data residual')
            xlabel(sprintf('Data #%d',id))
            ylabel('residual/noise realization')
            %title('')
            ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            print_mul(sprintf('%s_id%d_hist',fname,id))
        end
        
    end
    
catch
    keyboard
    cd(cwd);
    %close(f_handle)
    fprintf('%s : Cannot plot data response. \n',mfilename)
end


%% GO BACK TO STARTING DIRECTORY
cd(cwd)
