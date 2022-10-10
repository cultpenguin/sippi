% sippi_plot_loglikelihood Plot loglikelihood time series
%
% Call :
%    acc=sippi_plot_loglikelihood(logL,i_acc,N,itext)
%
function [acc,p,p2]=sippi_plot_loglikelihood(logL,i_acc,N,itext);

p2=[];
acc=NaN;
%cla;

nit=length(logL);
i=1:length(logL);

if nargin<2
    i_acc=1:1:length(logL);
end

if nargin<3
    N=1;
end

if nargin<4
    itext=max(i_acc)/5;
end

xlim=[1 nit];
p=plot(i,logL,'k-');
for ip=1:length(p);
    if ip==1; set(p(ip),'Color',[0 0 0]);end
    if ip==2; set(p(ip),'Color',[1 0 0]);end
    if ip==3; set(p(ip),'Color',[0 0 1]);end
    if ip==4; set(p(ip),'Color',[0 1 0]);end
    if ip==5; set(p(ip),'Color',[1 0 1]);end
end
%semilogy(i_acc,logL,'k-')
set(gca,'xlim',xlim)
grid on

%% logL=-N/2 +- sqrt(N/2)
if nargin>2   
    hold on
    p2=plot(xlim,[-1 -1].*N/2,'r-','linewidth',2);
    plot(xlim,[-1 -1].*N/2+2*sqrt(N/2),'r--');
    plot(xlim,[-1 -1].*N/2-2*sqrt(N/2),'r--');
    hold off
    % SCALE Y_AXIS
    n_int=ceil(0.01*length(logL));
    i1=(max([1,length(logL)-n_int]));
    y_min=min(logL(i1:end));
    try
        y_min=y_min-sqrt(N/2);
    end
    ylim=get(gca,'ylim');
    %if (y_min(1)*2)<ylim(1)
        ylim(1)=y_min(1)*2;
    %end
    set(gca,'ylim',ylim);
else
    return
end



%% Perhaps add text that indicates acceptance rate ?

% ylim=get(gca,'ylim');
% dy=diff(ylim);
% ylim(1)=ylim(1)-.1*dy;
% ylim(2)=ylim(2)+.1*dy;
% set(gca,'ylim',ylim);
% ntext=min([ceil(length(logL)/30) 15]);
% ntext=ceil(length(logL)/itext);
% ii=unique(round(linspace(1,length(logL),ntext)));
% dt=.1.*(ylim(2)-ylim(1));
% if length(ii)>2;
%     for i=2:(length(ii))
%         j=ii(i);
%         acc=(ii(i)-ii(i-1)) / (i_acc(ii(i))-i_acc(ii(i-1)));
%         y=min([logL(j)+dt,ylim(2)]);
%         %t=text(i_acc(j),min([logL(j)+dt,ylim(2)]),sprintf('%2.1f',100.*acc));
%         %set(t,'FontSize',8,'HorizontalAlignment','center')
%     end
% end

%set(gca,'xlim',xlim)
%set(gca,'ylim',ylim)
