% sippi_plot_current_model;
%
% Call :
%   sippi_plot_current_model(mcmc,data,d,m_current,prior);
%
function sippi_plot_current_model(mcmc,data,d,m_current,prior);

col=[
    0 0 0 
    1 0 0
    0 1 0 
    0 0 1
    1 1 0
    0 0 1
    .5 .5 .5
    ];


%% PLOT CURRENT MODELS
if nargin>3
    nm=length(m_current);
    figure_focus(100);subplot(1,1,1);
    sippi_plot_model(prior,m_current,1:nm,1,100);
    
    if isfield(mcmc,'m_ref');
        figure_focus(110);subplot(1,1,1);
        sippi_plot_model(prior,mcmc.m_ref,1:nm,1,110);
    end    
      
end

%% PLOT DATA RESPONSE
try
    sippi_plot_data(d,data);
end

%%

figure_focus(3);
%%
subplot(1,3,1);
try
    % CHANGE TO SUPPORT MORE THAN ONE DATA SET
    sippi_plot_loglikelihood(mcmc.logL(1:mcmc.i),mcmc.acc(1:mcmc.i),data{1}.N);
end
if exist('logL_ref','var');
    xlim=get(gca,'xlim');
    hold on
    plot(xlim,[1 1].*logL_ref,'k--','linewidth',2)
    hold off
end
%%
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
% CHANGE TO SUPPORT MORE THAN ONE DATA SET
%hx=linspace(-5,5,60);

for id=1:length(data);
    %try
        if isfield(data{id},'Cd')
            [h,hx]=hist(data{id}.d_obs(data{id}.i_use)-d{id}(:),30);
            h_true=normpdf(hx,0,sqrt(data{id}.Cd(1)));
        elseif isfield(data{id},'d_std')
            [h,hx]=hist(data{id}.d_obs(data{id}.i_use)-d{id}(:),30);
            h_true=normpdf(hx,0,sqrt(data{id}.d_std(1).^2));
        else
            [h,hx]=hist(data{id}.d_obs(data{id}.i_use)-d{id}(:),30);
            h_true=normpdf(hx,0,sqrt(data{id}.d_var(1)));
        end
    %end
    %h_true=normpdf(hx,0,sqrt(data{id}.Cd(1)));
    h_true=sum(h)*h_true./sum(h_true);
    bar(hx,[h],'EdgeColor',col(id+1,:),'FaceColor','none','BarWidth',1)
    if (id==1);hold on;end
    plot(hx,[h],'Color',col(id+1,:),'LineWidth',2)
    plot(hx,[h_true],'color',col(id+1,:).*.8,'LineWidth',2)
end
hold off
xlabel('d_{obs}-d_{cur}')
set(gca,'xlim',[min(hx) max(hx)])
%set(gca,'ylim',[0 1.2*max(h_true)])
drawnow;