function sippi_plot_prior(prior,m,im_array,use_colorbar,fhandle);
% sippi_plot_prior Plot a 'model', i.e. a realization of the prior model
%
%
% Call :
%   sippi_plot_prior(prior,m,im_array);
%
%   prior : Matlab structure for SIPPI prior model
%   m : Matlab structure for SIPPI realization
%   im_array : integer array of type of models to plot (typically 1)
%
%
% Example
%   m=sippi_prior(prior);
%   sippi_plot_prior(prior,m);
%
%   m=sippi_prior(prior);
%   sippi_plot_prior(prior,m,2);
%
% See also sippi_plot_prior, sippi_prior

if ischar(prior);
    disp(sprintf('%s: You probably intend to call sippi_plot_prior_sample',mfilename));   
    sippi_plot_prior_sample(prior);
    return
end

if nargin>1
    if ~iscell(m);
    disp(sprintf('%s: You probably intend to call sippi_plot_prior_sample',mfilename));   
    sippi_plot_prior_sample(prior);
    return
    end
end

% Check for initialization
for im=1:length(prior);
    if ~isfield(prior{im},'init')
        prior=sippi_prior_init(prior);
    end
end

if nargin<2 m=sippi_prior(prior); end
if nargin<3, im_array=1:1:length(m);end
if nargin<4, use_colorbar=1;end
if nargin<5, f_id=1;end
if nargin<3, subp=[1 1 1];end

if ~iscell(m)
    n{1}=m;
    m=n;
    clear n;
end

for im=im_array;

    dim=prior{im}.dim;
    x=prior{im}.x;
    y=prior{im}.y;
    z=prior{im}.z;
    
    nxyz=find(dim>1);
    ndim=length(nxyz);
    
    
    if ndim>0
        
        if isfield(prior{im},'cax');
            cax=prior{im}.cax;
        else
            cax=[min(m{im}(:)) max(m{im}(:))];
        end
        
        %f1=figure_focus;
        %figure_focus(f1+im-1);
        if nargin<5
            figure_focus(99+im);
        else
            try
            figure_focus(fhandle);
            end
        end
        
        if ndim==3;
            %% 3D MAP
            xslice = [max(x) (min(x)+max(x))/2];
            yslice = max(y);
            zslice = [ (min(z)+max(z))/2 max(z)];
            s=slice(prior{im}.x,prior{im}.y,prior{im}.z,m{im},xslice,yslice,zslice);
            for i=1:length(s)
                set(s(i),'Edgealpha',0);
            end
            %alpha(.85)
            view([-40 30])
            axis image
            xlabel('X');ylabel('Y');zlabel('Z')
            set(gca,'xdir','normal')
            set(gca,'ydir','normal')
            set(gca,'zdir','revers')
            doPlotIso=0;
            if doPlotIso==1
                p=patch(isosurface(x,y,z,m{im},cax(1)+.5*diff(cax)));
                isonormals(x,y,z,m{im},p);
                set(p,'FaceColor','red','EdgeColor','none');
                daspect([1 1 1])
                view(3); axis tight
                camlight
                camlight(160,10)
                %camlight('right')
                %camlight('left')
                lighting gouraud
                %lighting phong
                shading flat
                
            end
            
            if isfield(prior{im},'daspect');
                daspect(prior{im}.daspect);
            end
            
            try;caxis(cax);end
            
            colormap(sippi_colormap);
            if use_colorbar==1
                %colorbar_shift;
                colorbar;
            end
            
            
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        
        elseif ndim==2
            imagesc(x,y,squeeze(m{im}));
            %shading interp;
            axis image
            
            if isfield(prior{im},'daspect');
                daspect(prior{im}.daspect);
            end
            try;caxis(cax);end
            
            colormap(sippi_colormap);
            if use_colorbar==1
                colorbar_shift;
            end
            xlabel('X');
            ylabel('Y');
            if isfield(prior{im},'ydir');
                set(gca,'ydir',prior{im}.ydir);
            end            
            
        elseif ndim==1
            try
                if length(x)>length(y);
                    try
                        plot(x,m{im},'k-*')
                        xlabel('X');
                    catch
                        plot(m{im},'k-*')
                        xlabel('X [#]');                        
                    end
                else
                    plot(y,m{im},'k-*')
                    xlabel('Y');
                end
                ylabel(prior{im}.name)
                try;set(gca,'ylim',cax);end
            catch
                sippi_verbose(sprintf('%s could not plot model #%d',mfilename,im),1)
            end
        end
        
        if isfield(prior{im},'name');
            title(prior{im}.name,'interpreter','none')
        end
    else
        try;sippi_verbose(sprintf('%25s : m_%d(1)=%g',prior{im}.name,im,m{im}(1)));end
    end
end