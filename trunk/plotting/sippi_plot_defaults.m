% sippi_plot_defaults: Sets default options for ploting (such as fontsize)
%
% Call : 
%   options==sippi_plot_defaults(options);
%
%   % ALWAYS USE DEFULT SETTING (overrules options.axis) 
%   overrule=1; % {default overrule=0)
%   options==sippi_plot_defaults(options,overrule);
%   
%
function options=sippi_plot_defaults(options,overrule);
if nargin<2
    overrule=0;
end
if overrule==1;
    if isfield(options,'axis');
        options=rmfield(options,'axis');    
    end
end

options.axis.null='';
if ~isfield(options.axis,'fontsize');options.axis.fontsize=20;end
if ~isfield(options.axis,'width');options.axis.width=8;end
if ~isfield(options.axis,'height');options.axis.height=8;end
if ~isfield(options.axis,'w0');options.axis.w0=2;end
if ~isfield(options.axis,'h0');options.axis.h0=2;end
