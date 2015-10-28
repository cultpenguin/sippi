% sippi_plot_defaults: Sets default options for ploting (such as fontsize)
%
% Call :
%   options==sippi_plot_defaults(options);
%
%   % ALWAYS USE DEFULT SETTING (overrules options.axis)
%   overrule=1; % {default overrule=0)
%   options==sippi_plot_defaults(options,overrule);
%
% See also: sippi_plot_posterior, sippi_plot_posterior_2d_marg
%
function options=sippi_plot_defaults(options,overrule);
if nargin<2
    overrule=0;
end
if overrule==1;
    if isfield(options,'plot');
        options=rmfield(options,'plot');
    end
end

options.plot.null='';

options.plot.axis.null='';

%% AXES
options.axis.null='';
if ~isfield(options.plot.axis,'fontsize');options.plot.axis.fontsize=20;end
if ~isfield(options.plot.axis,'width');options.plot.axis.width=8;end
if ~isfield(options.plot.axis,'height');options.plot.axis.height=8;end
if ~isfield(options.plot.axis,'w0');options.plot.axis.w0=2;end
if ~isfield(options.plot.axis,'h0');options.plot.axis.h0=2;end

%% 'POSTERIOR' DATA 
options.plot.data.null='';
if ~isfield(options.plot.data,'show_max')
    % options.plot.data.show_max=Inf;
    options.plot.data.show_max=1;
end


%% 2D MARG
options.plot.marg2d.null='';
if ~isfield(options.plot.marg2d,'pl_marg2d_scatter')
    options.plot.marg2d.pl_marg2d_scatter=0;
end
if ~isfield(options.plot.marg2d,'pl_marg2d_image')
    options.plot.marg2d.pl_marg2d_image=0;
end
if ~isfield(options.plot.marg2d,'pl_marg2d_scatter_combine')
    options.plot.marg2d.pl_marg2d_scatter_combined=1;
end
if ~isfield(options.plot.marg2d,'pl_marg2d_image_combined')
    options.plot.marg2d.pl_marg2d_image_combined=1;
end
% use HPD for image plots
if ~isfield(options.plot.marg2d,'pl_marg2d_hpd')
%    options.plot.marg2d.pl_marg2d_hpd=0; % use scatter type instead
    options.plot.marg2d.pl_marg2d_hpd=1; % use HPD
end

if ~isfield(options.plot.marg2d,'NX'),options.plot.marg2d.NX=21;end
if ~isfield(options.plot.marg2d,'NY'),options.plot.marg2d.NY=21;end

if ~isfield(options.plot.marg2d,'hpd_interval'),options.plot.marg2d.hpd_interval=[.01,.1,.5,.9];end

%% use suptitle?
if ~isfield(options.plot,'suptitle'),options.plot.suptitle=0;end

%% SAMPLE TYPE. Skip or keep realization obatined using seqeuntial Gibbs sample
if ~isfield(options.plot,'skip_seq_gibbs'),
    options.plot.skip_seq_gibbs=1;
end


%% color codes
if ~isfield(options.plot,'color_codes');
    options.plot.color_codes=[
        0 0 0
        1 0 0
        0 1 0%% 2D MARG

        0 0 1
        1 1 0
        0 0 1
        .5 .5 .5
        ];
end

%% HARDCOPY types using print_mul, see print_mul for examples
if ~isfield(options.plot,'hardcopy_types');
    %options.plot.hardcopy_types=0; % no hardcopy
    options.plot.hardcopy_types=1; % PNG hardcopy
    %options.plot.hardcopy_types=2; % PDF hardcopy
end

%% SHOW_MAXIMUM NUMBER OF DATA IN SIPPI_PLOT_DATA
if ~isfield(options.plot,'plot_data_max_data');
    options.plot.plot_data_max_data=5;
end

