function sippi_plot_movie(fname,im_array,n_frames,skip_burnin,i_chain);
% sippi_plot_movie plot movie of prior and posterior realizations
%
% Call :
%   sippi_plot_movie(fname);
%   sippi_plot_movie(fname,im_array,n_frames,skip_burnin,i_chain);
%      fname : name of folder with results (e.g. options.txt)
%      im_array : array of indexes of model parameters to make into movies
%      n_frames [200] : number of frames in movie
%      skip_burnin [200] : start movie after burn_in;
%      i_chain[1]: make movie of chain number 'i_chain' (new 22/05/2014)
%
% Ex:
% sippi_plot_movie('20130812_Metropolis');
% sippi_plot_movie(options.txt);
%
% %% 1000 realization including burn-in, for prior number 1
% sippi_plot_movie('20130812_Metropolis',1,1000,0);
%
% Using options.plot.skip_seq_gibbs=1, (set in sippi_plot_defaults)
% removes realizations obtained using sequential Gibbs sampling 
% (equivalent to setting skip_burnin=1)
%
% See also: sippi_plot_defaults
%

% DOES NOT WORK WITH SIPPI_LEAST_SQUARES
% TEST THAT IT WORKS WITH sippi_rejection

if nargin==0;
    [f1,fname]=fileparts(pwd);
end
if nargin<3, n_frames=400;end

options=sippi_plot_defaults;

if nargin<4,
    skip_burnin=options.plot.skip_seq_gibbs;
end

if nargin<5,
    i_chain=1; % make movie of chain 1 realizations
end
n_frames_org=n_frames;

cwd=pwd;

%% DATA
if ischar(fname)
    try
        cd(fname);
        % load([fname,'.mat'],'prior','options','mcmc');
        load([fname,'.mat'],'prior','mcmc');
    catch
        % load([fname,'.mat'],'prior','options','mcmc');
        load([fname,'.mat'],'prior','mcmc');
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
if ~isfield(options,'txt')
    [p,options.txt]=fileparts(plotdir);
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

FrameRate = 30;
Quality = 95;
%%
for im=im_array
    
    ndim=length(find(prior{im}.dim>1));
    if ndim>0 % ONLY PLOT MOVIE FOR 1D, 2D and 3D PARAMETERS
        
        if ~isfield(mcmc,'i_sample');
            mcmc.i_sample=1;
        end
        
        N=prod(prior{im}.dim);
        if skip_burnin
            i1=ceil(prior{im}.seq_gibbs.i_update_step_max/mcmc.i_sample);
        else
            i1=1;
        end
        
        ns=mcmc.nite/mcmc.i_sample;
            
        n_frames=min([n_frames (ns-i1+1)]);
        i_frames=ceil(linspace(i1,ns,n_frames));
        if n_frames>0;
            
            try;disp(sprintf('Prior #%d : %s, plotting %d of %d posterior realizations',im,prior{im}.name,n_frames,ns));end
                        
            %% POSTERIOR
            vname=sprintf('%s_m%d_posterior_C%d',options.txt,im,i_chain);
            try
                if exist(vname,'file');
                    delete(vname);
                end
            end
            
            FrameRate = 25;
            Quality = 95;
            
            writerObj = VideoWriter(vname);
            %writerObj = VideoWriter(vname,'MPEG-4'); % Awful quality ?
            writerObj.FrameRate=FrameRate;
            writerObj.Quality=Quality;
            open(writerObj);
            
            
            fname=sprintf('%s_m%d_C%d.asc',options.txt,im,i_chain);
            if ~exist(fname,'file')
                % non-chains formatted output
                fname=sprintf('%s_m%d.asc',options.txt,im);
            end
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
                    
                    sippi_plot_prior(prior,m,im);
                    text(.02,.02,sprintf('#%05d, posterior',i),'units','normalized')
                    drawnow;
                    frame = getframe(gcf);
                    try
                        writeVideo(writerObj,frame);
                    catch
                        disp(sprintf('%s: problems writing %s',mfilename,fname));
                    end
                end
            end
            close(writerObj);
            fclose(fid);
            
        else
            try;disp(sprintf('Prior #%d : %s, no ''posterior'' frames to plot.',im,prior{im}.name));end
            i_frames=1:1:n_frames_org;
            
        end
        
        %% PRIOR
        prior{im}.seq_gibbs.step=prior{im}.seq_gibbs.step_max;
        prior{im}.perturb=0;
        vname=sprintf('%s_m%d_prior_C%d',options.txt,im,i_chain);
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
        
        for i=1:length(i_frames)
            [m]=sippi_prior(prior);
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
