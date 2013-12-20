function sippi_plot_movie(fname,im_array,n_frames,skip_burnin);
% sippi_plot_movie plot movie of prior and posterior realizations
%
% Call :
%   sippi_plot_movie(fname);
%   sippi_plot_movie(fname,im_array,n_frames,skip_burnin);
%      fname : name of folder with results (e.g. options.txt)
%      im_array : array of indexes of model parameters to make into movies
%      n_frames [200] : number of frames in movie
%      skip_burnin [200] : start movie after burn_in;
%
% Ex: 
% sippi_plot_movie('20130812_Metropolis');
% sippi_plot_movie(options.txt);
%
% %% 1000 realization including burn-in, for prior number 1
% sippi_plot_movie('20130812_Metropolis',1,1000,0);
%


if nargin==0;
    [f1,fname]=fileparts(pwd);
end

if nargin<4, skip_burnin=1;end
if nargin<3, n_frames=100;end


cwd=pwd;

%% DATA
if isstr(fname)
    try
        cd(fname);
        load([fname,'.mat']);
    catch
        load([fname,'.mat']);
    end
    
else
    disp(':/')
    return
end

plotdir=pwd;
try
    fname=options.txt;
end

if ~isfield(options,'FS')
    options.FS=12;
end

if nargin<2
    im_array=1:1:length(prior);
end

%%
options.axis_fontsize=8;
options.width=10;
options.height=10;
options.w0=2;
options.h0=2;

%%
for im=im_array
    
    ndim=length(find(prior{im}.dim>1));
    if ndim>1 % ONLY PLOT MOVIE FOR 2D and 3D PARAMETERS
        
        N=prod(prior{im}.dim);
        
        if skip_burnin
            i1=ceil(prior{1}.seq_gibbs.i_update_step_max/mcmc.i_sample);
        else
            i1=1;
        end
        ns=mcmc.nite/mcmc.i_sample;
        n_frames=min([n_frames (ns-i1)]);
        i_frames=ceil(linspace(i1,ns,n_frames));
        

        try;disp(sprintf('Prior #%d : %s, plotting %d of %d posterior realizations',im,prior{im}.name,n_frames,ns));end

        
        %% POSTERIOR
        vname=sprintf('%s_m%d_posterior',options.txt,im);
        try
            if exist(vname,'file');
                delete(vname);
            end
        end
        
        writerObj = VideoWriter(vname);
        %writerObj = VideoWriter(vname,'MPEG-4'); % Awful quality ?
        writerObj.FrameRate=25;
        open(writerObj);
        
        fname=sprintf('%s_m%d.asc',options.txt,im);
        fid=fopen(fname,'r');
        
        i=0;
        while ~feof(fid);
            i=i+1;
            d=fscanf(fid,'%g',N);
            
            if ~isempty(find(i==i_frames))
                if prior{im}.dim(3)>1
                    % 3D
                    real=reshape(d,length(prior{im}.y),length(prior{im}.x),length(prior{im}.z));
                elseif prior{im}.dim(2)>1
                    % 2D
                    real=reshape(d,length(prior{im}.y),length(prior{im}.x));
                else
                    % 1D
                    real=d;
                end
                m{im}=real;
                
                sippi_plot_model(prior,m,im);
                text(.02,.02,sprintf('#%05d',i),'units','normalized')
                drawnow;
                frame = getframe;
                writeVideo(writerObj,frame);
            end
        end
        close(writerObj);
        fclose(fid);
        
        %% PRIOR
        for i=1:length(prior);
            prior{i}.seq_gibbs.step=prior{i}.seq_gibbs.step_max;
        end
        vname=sprintf('%s_m%d_prior.mp4',options.txt,im);
        try
            if exist(vname,'file');
                delete(vname);
            end
        end
        
        writerObj = VideoWriter(vname);
        %writerObj = VideoWriter(vname,'MPEG-4'); % Awful quality ?
        writerObj.FrameRate=25;
        writerObj.Quality=100;
        open(writerObj);
        
        for i=1:length(i_frames)
            m=sippi_prior(prior,m);
            sippi_plot_model(prior,m,im);
            text(.02,.02,sprintf('#%05d',i),'units','normalized')
            drawnow;
            frame = getframe;
            writeVideo(writerObj,frame);
        end
        close(writerObj);
    
        
    else
        try;disp(sprintf('''%s'' only applies to 2D/3D parameters (not prior #%d : %s)',mfilename,im,prior{im}.name));end
    end
end
%%
cd(cwd)
