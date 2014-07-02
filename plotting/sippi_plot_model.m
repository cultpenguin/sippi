function sippi_plot_model(prior,m,im_array,use_colorbar,fhandle);

disp(sprintf('''%s'' is obsolete, please use ''%s'' instead',mfilename,'sippi_plot_prior'))
if nargin==1
    sippi_plot_prior(prior);
elseif nargin==2
    sippi_plot_prior(prior,m)
elseif nargin==3
    sippi_plot_prior(prior,m,im_array)
elseif nargin==4
    sippi_plot_prior(prior,m,im_array,use_colorbar);
elseif nargin==5
    sippi_plot_prior(prior,m,im_array,use_colorbar,fhandle);
end
