% sippi_prior_snesim : SNESIM type Gaussian prior for SIPPI
%
%                      Using SNESIM form
%                      https://github.com/SCRFpublic/snesim-standalone
%                      Please remember to recompile SNESIM to uou needs, 
%                      before using it with SIPPI
%
%
%% Example:
%    ip=1;
%    prior{ip}.type='snesim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.ti=channels;
%    % prior{ip}.ti=maze;
%
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%    figure(1);imagesc(prior{ip}.ti);axis image
%
%% Example: scaling and rotation
%    ip=1;
%    prior{ip}.type='snesim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.ti=channels;
%    prior{ip}.scaling=[.1];
%    prior{ip}.rotation=[10];
%
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%    figure(1);imagesc(prior{ip}.ti);axis image
%
%% Hard data
%   % hard data are given using either matrix of 4 columns (x,y,z,val)
%   % or as a 4 column EAS file (x,y,z,val)
%   d_hard=[1 1 0 0; 1 2 0 0; 2 2 0 1 ];
%   prior{ip}.hard_data=d_hard;
%
%   write_eas('snesim_hard.dat',d_hard);
%   prior{ip}.hard_data='snesim_hard.dat';
%
%% Soft data
%   % soft mush be provided as a matrix of the same size as the simulation
%   % grid
%   d_soft(:,:,1)=ones(80,80).*NaN
%   d_soft(:,:,2)=1-d_soft(:,:,1);
%   prior{ip}.soft_data_grid=d_soft;
%
%
%   % Optionally the soft data can be provided as point data, in which case
%   % a grid, the size of the simulation grid, with soft data values will be computed
%   d_soft=[1 1 0 0.2 0.8; 1 2 0 0.1 0.9; 2 2 0  0.05 0.95];
%   prior{ip}.soft_data=d_soft;
%
%% Sequential Gibbs sampling type 1 (box selection of pixels)
%    prior{ip}.seq_gibbs.type=1;%
%    prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
%    [m,prior]=sippi_prior(prior);
%    for i=1:10;
%       [m,prior]=sippi_prior(prior,m);
%       sippi_plot_prior(prior,m);
%       drawnow;
%    end
%
%% Sequential Gibbs sampling type 2 (random pixels)
%    prior{ip}.seq_gibbs.type=2;%
%    prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
%    [m,prior]=sippi_prior(prior);
%    for i=1:10;
%       [m,prior]=sippi_prior(prior,m);
%       sippi_plot_prior(prior,m);
%       drawnow;
%    end
%
% See also: sippi_prior, ti
%
function [m_propose,prior]=sippi_prior_snesim(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end


% SNESIM

%% REMOVE CONDITIONAL DATA.
%% FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
%if isfield(prior{ip}.S,'f_obs')
%    prior{ip}.S=rmfield(prior{ip}.S,'f_obs');
%end
%prior{ip}.S.XML.parameters.Hard_Data.grid='';
%prior{ip}.S.XML.parameters.Hard_Data.property='';

% force nsim=1 in sippi
prior{ip}.S.nsim=1;

% set random seed
if isfield(prior{ip},'seed');
    if (prior{ip}.seed)==0
        prior{ip}.S.rseed=ceil(rand(1).*1e+6);
    else
        prior{ip}.S.rseed=prior{ip}.seed;
    end
else
    prior{ip}.S.rseed=ceil(rand(1).*1e+6);
end

% optionally set rotation and affinity
set_aff=0;
if isfield(prior{ip},'rotation')|isfield(prior{ip},'scaling')
    set_aff=1;
end
if set_aff==1
    if ~isfield(prior{ip},'rotation'), prior{ip}.rotation=1; end
    if ~isfield(prior{ip},'scaling'), prior{ip}.scaling=1; end
    prior{ip}.S=snesim_set_rotation_affinity(prior{ip}.S,prior{ip}.rotation,1./prior{ip}.scaling);
    prior{ip}.S.frotaff.use=1;
end

%% optionally set hard data
if isfield(prior{ip},'hard_data');
    if ischar(prior{ip}.hard_data)
        % Hard data is provided in file
        prior{ip}.S.fconddata.fname=prior{ip}.hard_data;
    else
        % save hard data, and set hard data filename
        prior{ip}.S.fconddata.fname='snesim_hard.dat';
        filename=prior{ip}.S.fconddata.fname;
        sippi_verbose(sprintf('%s: saving hard data to %s',mfilename,filename));
        write_eas(filename,prior{ip}.hard_data);
    end
else
    if exist('snesim_hard_dummy.dat');
        try;delete('snesim_hard_dummy.dat');end;
    end
    prior{ip}.S.fconddata.fname='snesim_hard_dummy.dat';
end

%% optionally set soft data
if isfield(prior{ip},'soft_data');
    if ischar(prior{ip}.soft_data)
        % soft data is provided in file
        prior{ip}.S.flocalprob.fname=prior{ip}.soft_data;
    else
        if ~isfield(prior{ip},'soft_data_grid');
            % only update to grid if grid does not exist
            sippi_verbose(sprintf('%s: converting soft data points to grid',mfilename));
            % compute soft data grid from point data
            ncat=length(unique(prior{ip}.ti(:)));
            for ic=1:ncat
                if prior{ip}.ndim==1
                    soft_data_grid(:,ic)=(1/ncat)+prior{1}.xx.*0;
                elseif prior{ip}.ndim==2
                    soft_data_grid(:,:,ic)=(1/ncat)+prior{1}.xx.*0;
                else
                    soft_data_grid(:,:,:,ic)=(1/ncat)+prior{1}.xx.*0;
                end
            end
            % mv point data into grid
            for i=1:size(prior{ip}.soft_data,1);
                ix=find(prior{ip}.x==prior{ip}.soft_data(i,1));
                iy=find(prior{ip}.y==prior{ip}.soft_data(i,2));
                iz=find(prior{ip}.z==prior{ip}.soft_data(i,3));
                for ic=1:ncat
                    if prior{ip}.ndim==1
                        soft_data_grid(ix)=prior{ip}.soft_data(i,3+ic);
                    elseif prior{ip}.ndim==2
                        soft_data_grid(iy,ix,ic)=prior{ip}.soft_data(i,3+ic);
                    else
                        soft_data_grid(ix,iy,iz,ic)=prior{ip}.soft_data(i,3+ic)
                    end
                end
            end
            
            prior{ip}.soft_data_grid=soft_data_grid;
            
        end
    end
end

%% optionally set soft data grid
if isfield(prior{ip},'soft_data_grid');
    if ischar(prior{ip}.soft_data_grid)
        % soft data is provided in file
        prior{ip}.S.flocalprob.fname=prior{ip}.soft_data_grid;
    else
        % save soft data, and set hard data filename
        filename=prior{ip}.S.flocalprob.fname;
        sippi_verbose(sprintf('%s: saving soft data to %s',mfilename,filename));
        if prior{ip}.ndim==1;
            write_eas(filename,prior{ip}.soft_data_grid);
        elseif prior{ip}.ndim==2;
            ncat=length(unique(prior{ip}.ti(:)));
            for i=1:(ncat);
                p=prior{ip}.soft_data_grid(:,:,i)';
                soft_data(:,i)=p(:);
            end
            write_eas(filename,soft_data);
        else
            sippi_verbose(sprintf('%s: soft_data_grid not implemented for 3D data',mfilename));
        end
    end
    prior{ip}.S.condition_to_lp=1;
    prior{ip}.S.iauto=1;
    
else
    if exist('snesim_soft_dummy.dat');
        try;delete('snesim_soft_dummy.dat');end;
    end
    prior{ip}.S.flocalprob.fname='snesim_soft_dummy.dat';
    prior{ip}.S.condition_to_lp=0;
    prior{ip}.S.iauto=0;
end



%% sequential Gibbs resampling
if nargin>1
    
    % Convert values to indexes
    if isfield(prior{ip},'index_values');
        m = zeros(size(m_current{ip}))-1;
        for i=1:length(prior{ip}.index_values)
            try
                m(find(m_current{ip}==prior{ip}.m_values(i)))=prior{ip}.index_values(i);
            end
        end
        m_current{ip}=m;
    end
    
    % SEQUENTIAL GIBBS    
    sippi_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
    %prior{ip}.S=snesim_set_resim_data(prior{ip}.S,prior{ip}.S.D,[10 10]);    
    prior{ip}.S=snesim_set_resim_data(prior{ip}.S,m_current{ip},prior{ip}.seq_gibbs.step,prior{ip}.seq_gibbs.type);
    
end

%% RUN SNESIM
prior{ip}.S = snesim(prior{ip}.S,prior{ip}.x,prior{ip}.y,prior{ip}.z);
m_propose{ip} = prior{ip}.S.D;

%% Convert indexes to values
if isfield(prior{ip},'index_values');
    for i=1:length(prior{ip}.index_values)
        m_propose{ip}(find(m_propose{ip}==prior{ip}.index_values(i)))=prior{ip}.m_values(i);
    end
end
