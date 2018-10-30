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
% A maximum of options.plot.plot_data_max_data (def:=5) data is plotted.
% Change the value in sippi_plot_defaults.m
%
% A specific m-file for handling plotting for a specfic type of data can
% be implemented, and can be called instead of sippi_plot_data, if the name
% of the m-file is set in the options.sippi_plot_data_functions is set, as
% e.g:
%    options.sippi_plot_data_function='sippi_plot_data_gpr';
% then 
%    sippi_plot_data(d,data,[],options)
%    sippi_plot_data(d,data,1:length(d),options); 
% will be equivalent to call
%    sippi_plot_data_gpr(d,data); 
%
% 
%
% 
% See also sippi_plot_defaults, sippi_plot_data_gpr
%
function sippi_plot_data(d,data,id_arr,options);

if nargin<3;
  id_arr=1:length(d);
end
if isempty(id_arr);
  id_arr=1:length(d);
end    

options.null='';

if nargin>3
    if isfield(options,'sippi_plot_data_function');
        sippi_verbose(sprintf('%s : Using ''%s'' instead of sippi_plot_data',mfilename,options.sippi_plot_data_function),1);
        if nargin==1;
            feval(options.sippi_plot_data_function,d);
        elseif nargin==4            
            feval(options.sippi_plot_data_function,d,data,id_arr,options);
        else
            feval(options.sippi_plot_data_function,d,data,id_arr);
        end
        return
    end
end

options=sippi_plot_defaults(options);

if (length(id_arr)>options.plot.plot_data_max_data)
    sippi_verbose(sprintf('%s: Plotting only the first %d (of %d) data',mfilename,options.plot.plot_data_max_data,length(id_arr)));
    sippi_verbose(sprintf('%s: -> as defined in sippi_plot_defaults.m',mfilename));
    id_arr=id_arr(1:options.plot.plot_data_max_data);
end

for id=id_arr;
   
    figure_focus(30+id);
    
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
