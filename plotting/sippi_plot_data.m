% sippi_plot_data: Plot data response
%
% Call:
%    sippi_plot_data(d,data);
%
% sippi_plot_data provides a very simple way to plot data. 
% A more appropriate data plot can be implemented by implementing a new
% mfile called "sippi_plot_data" and add it the Matlab path before the 
% main SIPPI folders
%
function sippi_plot_data(d,data,id_arr);

if nargin<3;
  id_arr=1:lenghth(d);
end

for id=id_arr;
    
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
