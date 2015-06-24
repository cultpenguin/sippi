% sippi_plot_data: plot NMR data response
%                  overwrites sippi_plot_data plot data in SIPPI
%
% Call.
%    sippi_plot_data(d,data);
%
function sippi_plot_data(d,data);


for id=1:length(d);
    
    figure_focus(20+id);
    
    if ~isfield(data{id},'i_use')
        data{id}.i_use=1:1:length(data{id}.d_obs);
    end
    
    plot(d{id},'k-');
    %wiggle(d{id});
    if nargin>1
        hold on
        plot(data{id}.d_obs(data{id}.i_use),'r-');
        %wiggle(data{id}.d_obs);
        hold off
        legend('d_{cur}','d_{obs}')
    else
        legend('d_{obs}')
    end
    xlabel('data #')
    ylabel('')
end
