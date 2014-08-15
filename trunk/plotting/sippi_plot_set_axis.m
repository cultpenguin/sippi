% sippi_plot_set_axis
% see also sippi_plot_defaults

function sippi_plot_set_axis(options)
if nargin==1
    options=sippi_plot_defaults(options);
else
    options=sippi_plot_defaults;
end    

%ppp(options.plot.axis.width,options.plot.axis.height,options.plot.axis.fontsize,options.plot.axis.w0,options.plot.axis.h0);
            
set(gca,'FontSize',options.plot.axis.fontsize)

set(get(gca,'Xlabel'),'FontSize',options.plot.axis.fontsize+2)
set(get(gca,'Ylabel'),'FontSize',options.plot.axis.fontsize+2)
set(get(gca,'Zlabel'),'FontSize',options.plot.axis.fontsize+2)

