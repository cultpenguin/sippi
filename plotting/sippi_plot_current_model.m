% sippi_plot_current_model Plots the current model during Metropolis sampling
%
% Call :
%   sippi_plot_current_model(mcmc,data,d,m_current,prior,options);
%
function sippi_plot_current_model(mcmc,data,d,m_current,prior,options);
options.null='';


col=[
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    0 0 1
    .5 .5 .5
    .4 .4 .4
    .3 .3 .3
    .2 .2 .2
    ];


%% PLOT CURRENT MODELS
if nargin>3
    sippi_plot_prior(prior,m_current);
    
    % If reference model is set, then plot the reference model   
    if isfield(options.mcmc,'m_ref');
        try
        for ip=1:length(prior);
            sippi_plot_prior(prior,options.mcmc.m_ref,ip,1,20+ip-1); 
            title(sprintf('Reference Prior #%d',ip))
        end
        catch
            keyboard
        end
    end
    
end

%% PLOT DATA RESPONSE
try
    sippi_plot_data(d,data,1:length(data),options);
end

%%

figure_focus(3);
subfigure(2,2,2)
%%
subplot(1,3,1);
try
    
    % CHANGE TO SUPPORT MORE THAN ONE DATA SET
    N=0;
    for i=1:length(data);
        if isfield(data{i},'i_use');
            N=N+length(data{i}.i_use);
        else
            N=N+length(data{i}.d_obs);
        end
    end
    if ~isfield(mcmc,'i');
        if size(mcmc.acc,1)>1
            mcmc.i=max(find(sum(mcmc.acc)==1));
        else
            mcmc.i=max(find(mcmc.acc==1));
        end
    end
    sippi_plot_loglikelihood(mcmc.logL(1:mcmc.i),mcmc.acc(1:mcmc.i),N);
catch
    sippi_verbose(sprintf('%s: failed to plot log likelihood curve',mfilename))
end

subplot(1,3,2);
try
    % MAKE SURE THIS WORKS FOR MANY PRIOR MODELS
    %    plot(1:1:mcmc.i,mcmc.step(1:1:mcmc.i),'k-')
    %else
    semilogy(1:1:mcmc.i,mcmc.step(:,1:1:mcmc.i),'k-')
    %end
    ylabel('Step length')
    xlabel('Iteration number')
end
%%
subplot(1,3,3);

show_max_data=2; % NEEDS TO BE AN ADJUSTABLE PARAMETER
for id=1:[min([show_max_data length(data)])];
    %try
    try
        [h,hx]=hist(data{id}.d_obs(data{id}.i_use)-d{id}(:),30);
    catch
        [h,hx]=hist(data{id}.d_obs(data{id}.i_use)-d{id}(:)',30);
    end
    bar(hx,[h],'EdgeColor',col(id+1,:),'FaceColor','none','BarWidth',1)

    if isfield(data{id},'Cd')
        h_true=normpdf(hx,0,sqrt(data{id}.Cd(1)));
    elseif isfield(data{id},'d_std')
        h_true=normpdf(hx,0,sqrt(data{id}.d_std(1).^2));
    else
        h_true=normpdf(hx,0,sqrt(data{id}.d_var(1)));
    end
    h_true=sum(h)*h_true./sum(h_true);
    if (id==1);hold on;end
    plot(hx,[h],'Color',col(id+1,:),'LineWidth',2)
    plot(hx,[h_true],'color',col(id+1,:).*.8,'LineWidth',2)
end
hold off
xlabel('d_{obs}-d_{cur}')
set(gca,'xlim',[min(hx) max(hx)])
%set(gca,'ylim',[0 1.2*max(h_true)])
drawnow;