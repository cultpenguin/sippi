function sippi_plot_data_fdem1d(d,data,id_arr,options);


if nargin<3;
    id_arr=1:length(d);
end
if isempty(id_arr);
    id_arr=1:length(d); per
end

options.null='';

if length(d{1})==12;
    xarr=[395,        1822,        3262,        8199,       38760,      128760];
else
    xarr=1:6;
end

for id=id_arr;
    
    figure_focus(20+id);
    
    col1=[0 0 0];
    col2=[1 1 1].*.6;
    
    Nd=12; % number of data per sounding
    
    N=length(d{id});
    
    if N==Nd
        N2=N/2;
        d_real=d{id}(1:2:N);
        d_imag=d{id}(2:2:N);
        
        if nargin>1;
            d_real_obs=data{id}.d_obs(1:2:N);
            d_imag_obs=data{id}.d_obs(2:2:N);
            d_real_std=data{id}.d_std(1:2:N);
            d_imag_std=data{id}.d_std(2:2:N);
            
            p(1)=errorbar(xarr,d_real_obs,2*d_real_std,'k','LineWidth',2);
            hold on
            plot(xarr,d_real_obs,'k+','MarkerSize',12)
            
            p(2)=errorbar(xarr+0.01*xarr,d_imag_obs,2*d_imag_std,'color',col2,'LineWidth',2);
            plot(xarr,d_imag_obs,'+','MarkerSize',12,'color',col2)
            
        end
        
        plot(xarr,d_real,'k.','MarkerSize',38)
        hold on
        plot(xarr,d_real,'k-','linewidth',1)
        
        plot(xarr,d_imag,'o','color',col2,'MarkerSize',12,'LineWidth',3)
        plot(xarr,d_imag,'-','color',col2,'linewidth',1)
        hold off
        
        xlabel('Frequency (Hz) #')
        legend([p(1) p(2)],'In phase +- 2\sigma','Quadrature +- 2\sigma')
        set(gca,'XScale','log')
        xlim([100 500000])
    else
        % MULD
        % N~=Nd
        
        if nargin>1;
            D_obs=reshape(data{id}.d_obs,Nd,N/Nd);
            D_std=reshape(data{id}.d_std,Nd,N/Nd);
            e=errorbar(D_obs(1:2:12,:)',D_std(1:2:12,:)','k-');
            hold on
            e2=errorbar(D_obs(2:2:12,:)',D_std(2:2:12,:)','-','color',col2);
        end
        Nobs=size(D_obs,2);
        D=reshape(d{id},Nd,N/Nd);
        d_real=D(1:2:end,:);
        d_imag=D(2:2:end,:);
        xx=repmat(xarr(:),1,Nobs);
        xx=1:size(D,2);
        plot(xx,d_real,'k.','MarkerSize',18)
        hold on
        %plot(xx,d_real,'k-')
        %plot(xx,d_imag,'-','color',col2)
        plot(xx,d_imag,'o','color',col2,'MarkerSize',7,'LineWidth',2)
        hold off
        
        xlabel('Frequency (Hz) #')
        legend('In phase +- 2\sigma','Quadrature +- 2\sigma')
        
    end
    ylim([-200 4500])
    
end
grid on