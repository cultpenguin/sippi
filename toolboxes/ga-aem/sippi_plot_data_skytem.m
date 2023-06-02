% sippi_plot_data_gaaem: 
%                  
% Call.
%    sippi_plot_data_gaaem(d,data);
%    sippi_plot_data_gaaem(d,data,forward);
%
% See also: sipi_plot_data
%
function sippi_plot_data_gaaem(d,data,id_arr,forward);

if nargin<4, 
    forward.null=[];
end
if nargin<3, 
    id_arr=1;
end
if nargin<2,
    data=[];
end

for id=1;
    
    figure_focus(20+id);
    if isfield(forward,'LM')
        d1 = d{id}(1:forward.LM.nw);
        d2 = d{id}((forward.LM.nw+1):end);
        
        h1=plot(forward.LM.wt.centre,d1,'-r.','linewidth',1);
        hold on
        h2=plot(forward.HM.wt.centre,d2,'-b.','linewidth',1);
        hold off
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        %set(gca,'xlim',[1e-5  1e-2]);
        %set(gca,'ylim',[1e-15 1e-7]);
        xlabel('Time (s)');
        ylabel('Response (V/A.m^4)');
        box on;
        grid on
        legend([h1 h2],'LM Z','HM Z');
   
        if ~isempty(data)
            hold on
            d1_obs = data{id}.d_obs(1:forward.LM.nw);
            d2_obs = data{id}.d_obs((forward.LM.nw+1):end);
        
            h1_d=plot(forward.LM.wt.centre,d1_obs,'r*','linewidth',1);
            h2_d=plot(forward.HM.wt.centre,d2_obs,'b*','linewidth',1);
            
            
            try
            d1_std = data{id}.d_std(1:forward.LM.nw);
            d2_std = data{id}.d_std((forward.LM.nw+1):end);
            errorbar(forward.LM.wt.centre,d1_obs,2*d1_std,'color','r')
            errorbar(forward.HM.wt.centre,d2_obs,2*d2_std,'color','b')
            end
            
            hold off
        
        end
        
    else
        plot(d{id},'k-')
        set(gca,'yscale','log');
        set(gca,'ylim',[1e-15 1e-7]);
        xlabel('Data #');
        ylabel('Response (V/A.m^4)');
        if ~isempty(data)
            ii = 1:length(data{id}.d_obs);
            hold on
            errorbar(ii,data{id}.d_obs,data{id}.d_std,'color','r');
            hold off
        end
    end
    
    
end


