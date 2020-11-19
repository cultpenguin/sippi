% sippi_plot_prior_movie: creates a movie file a random walk in the prior
% 
% Call:
%    sippi_plot_prior_movie(prior,n_frames,options,im_array);
%
% See also: sippi_plot_movie
%
function sippi_plot_prior_movie(prior,n_frames,options,im_array);

if nargin<2
    n_frames=100;
end
if nargin<4
    im_array=1:length(prior);
end
options.null='';
if ~isfield(options,'txt');
    options.txt='prior';
end
options=sippi_plot_defaults(options);

cwd=pwd;

%%
options.axis_fontsize=8;
options.width=10;
options.height=10;
options.w0=2;
options.h0=2;

FrameRate=30;
Quality=90;
cwd=pwd;

[m,prior]=sippi_prior(prior);

%%
for im=im_array
    
    ndim=length(find(prior{im}.dim>1));
    if ndim>0 % ONLY PLOT MOVIE FOR 1D, 2D and 3D PARAMETERS
        
        
        %% PRIOR
        %prior{im}.seq_gibbs.step=prior{im}.seq_gibbs.step_max;
        prior{im}.perturb=1;
        vname=sprintf('%s_m%d',options.txt,im);
        try
            if exist(vname,'file');
                delete(vname);
            end
        end
        
        writerObj = VideoWriter(vname);
        %writerObj = VideoWriter(vname,'MPEG-4'); % Awful quality ?
        writerObj.FrameRate=FrameRate;
        writerObj.Quality=Quality;
        open(writerObj);
        
        for i=1:n_frames
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m,im);
            text(.02,.02,sprintf('#%05d, prior',i),'units','normalized')
            drawnow;
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        close(writerObj);
        
        
    else
        try;sippi_versboe(sprintf('''%s'' only applies to 2D/3D parameters (not prior #%d : %s)',mfilename,im,prior{im}.name),1);end
    end
end
%%db
cd(cwd)
