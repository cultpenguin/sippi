% sippi_prior_mps : prior based on MPS
%
%                      Using SNESIM/ENESIM FROM 
%                      https://github.com/ergosimulation/mpslib
% 
%% Example:
%    ip=1;
%    prior{ip}.type='mps';
%    prior{ip}.method='mps_snesim';
%    prior{ip}.x=1:1:80;
%    prior{ip}.y=1:1:80;
%    prior{ip}.ti=channels;
%    % prior{ip}.ti=maze;
%
%    m=sippi_prior(prior);
%    sippi_plot_prior(prior,m)
%    figure(1);imagesc(prior{ip}.ti);axis image
%
%
%%  The specific algorithm use for MPS simulation is defined in ytje 'method' field
%    prior{ip}.method='mps_snesim';         % default, same as 'mps_snesim_tree'
%    prior{ip}.method='mps_snesim_tree';
%    prior{ip}.method='mps_snesim_list';
%    prior{ip}.method='mps_genesim';
%
%    All properties for each algorithm are availale in the prior{ip}.MPS
%    field
%
% See also: sippi_prior, ti, mps_cpp, sippi_sequential_gibbs_resim

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
%
%
function [m_propose,prior]=sippi_prior_mps(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if ~isfield(prior{ip},'ti')
    prior{ip}.ti=channels;
end

if ~isfield(prior{ip},'method')
    prior{ip}.method='mps_snesim_tree';
    % prior{ip}.method='mps_snesim_list';
    % prior{ip}.method='mps_genesim';
end

if ~exist('mps_cpp.m','file')
    disp(sprintf('%s: =============================================================',mfilename))
    disp(sprintf('%s: You need to install MPSlib to use %s',mfilename,prior{ip}.method))
    disp(sprintf('%s: MPSlib --> https://github.com/ergosimulation/mpslib/releases',mfilename))    
    disp(sprintf('%s: =============================================================',mfilename))
end


prior{ip}.MPS.method=prior{ip}.method;
%prior{ip}.MPS.template_size=[9 9 1];
prior{ip}.MPS.nreal=1;

if ~isfield(prior{ip}.MPS,'parameter_filename')
    prior{ip}.MPS.parameter_filename=sprintf('%s_%04d.txt',prior{ip}.method,round(1000*rand(1)));;
end

if prior{ip}.ndim==1;
    SIM=NaN.*ones(prior{ip}.dim(1));    
elseif prior{ip}.ndim==1;
    SIM=NaN.*ones(prior{ip}.dim(2),prior{ip}.dim(1));    
else
    SIM=NaN.*ones(prior{ip}.dim(2),prior{ip}.dim(1),prior{ip}.dim(3));    
end

% initialize prior and set the x,y,z dimensions
prior=sippi_prior_init(prior,ip);
if ~isfield(prior{ip}.MPS,'simulation_grid_size')
    prior{ip}.MPS.simulation_grid_size=[length(prior{ip}.x) length(prior{ip}.y) length(prior{ip}.z)];
end
prior{ip}.MPS.origin=[prior{ip}.x(1) prior{ip}.y(1) prior{ip}.z(1)];
prior{ip}.MPS.grid_cell_size=[1 1 1];
if prior{ip}.dim(1)>1
    prior{ip}.MPS.grid_cell_size(1)=prior{ip}.x(2)-prior{ip}.x(1);
end
if prior{ip}.dim(2)>1
    prior{ip}.MPS.grid_cell_size(2)=prior{ip}.y(2)-prior{ip}.y(1);
end
if prior{ip}.dim(3)>1
    prior{ip}.MPS.grid_cell_size(3)=prior{ip}.z(2)-prior{ip}.z(1);
end

% random seed?
% set random seed
if isfield(prior{ip},'seed');
    prior{ip}.MPS.rseed=prior{ip}.seed;
else
    prior{ip}.MPS.rseed=ceil(rand(1).*1e+6);
end

%% Convert values to indexes
if nargin>1
    if isfield(prior{ip},'index_values');
        m = zeros(size(m_current{ip}))-1;
        for i=1:length(prior{ip}.index_values)
            try
                m(find(m_current{ip}==prior{ip}.m_values(i)))=prior{ip}.index_values(i);
            end
        end
        m_current{ip}=m;
    end
end    

%% Sequential gibbs resampling
if nargin>1
   d_cond=sippi_get_resim_data(m_current,prior,ip);    
   % set hard data
   prior{ip}.hard_data=d_cond;
end
   

%% hard data?
if nargin<2
    % for unconditional simulation force
    % removal of hard data
    if isfield(prior{ip},'hard_data');
        prior{ip}=rmfield(prior{ip},'hard_data');
    end
end
    
if isfield(prior{ip},'hard_data');
    if ischar(prior{ip}.hard_data)
        % Hard data is provided in file
        prior{ip}.MPS.hard_data_filename=prior{ip}.hard_data;
    else
        % save hard data, and set hard data filename
        if ~isfield(prior{ip}.MPS,'hard_data_filename')
            prior{ip}.MPS.hard_data_filename='f_cond.dat';
        end
        filename_hard=prior{ip}.MPS.hard_data_filename;
        if isempty(prior{ip}.hard_data)
            if exist([pwd,filesep,filename_hard],'file');
                delete(filename_hard);
            end
        else
            sippi_verbose(sprintf('%s: saving hard data to %s',mfilename,filename_hard),2);
            write_eas(filename_hard,prior{ip}.hard_data);
        end
    end
else
    f_hard_data_dummy='mps_hard_data_dummy.dat';
    if exist(f_hard_data_dummy,'file');
        try;           
            delete(f_hard_data_dummy);
        catch
            disp(sprintf('%s: Could not delete %s',mfilename,f_hard_data_dummy))
        end
    end
    prior{ip}.MPS.hard_data_filename=f_hard_data_dummy;
end

%% soft data?
if isfield(prior{ip},'soft_data');
    if ischar(prior{ip}.soft_data)
        % Hard data is provided in file
        prior{ip}.MPS.soft_data_filename=prior{ip}.soft_data;
    else
        % save soft data, and set  soft filename
        try
            filename_soft=prior{ip}.MPS.soft_data_filename;
        catch
            filename_soft='mps_soft.dat';
        end
        prior{ip}.MPS.soft_data_filename=filename_soft;
        sippi_verbose(sprintf('%s: saving soft data to %s',mfilename,filename_soft));
        write_eas(filename_soft,prior{ip}.soft_data);
    end
elseif isfield(prior{ip},'soft_data_grid');
    if ischar(prior{ip}.soft_data_grid)
        % Hard data is provided in file
        prior{ip}.MPS.soft_data_filename=prior{ip}.soft_data_grid;
    else
        % convert grid file into [x y z p0 p1 ...] file
        ncat=length(unique(prior{ip}.ti(:)));
        soft_data(:,1)=prior{ip}.xx(:);
        soft_data(:,2)=prior{ip}.yy(:);
        soft_data(:,3)=prior{ip}.zz(:);
        for ic=1:ncat
            if prior{ip}.ndim==1
                d=prior{ip}.soft_data_grid(:,ic);
            elseif prior{ip}.ndim==2
                d=prior{ip}.soft_data_grid(:,:,ic);
            else
                d=prior{ip}.soft_data_grid(:,:,:,ic);
            end
            soft_data(:,3+ic)=d(:);

        end
        filename_soft=prior{ip}.MPS.soft_data_filename;
        sippi_verbose(sprintf('%s: saving soft data grid to %s',mfilename,filename_soft));
        write_eas(filename_soft,soft_data);
        
     end
else
    if exist('mps_soft_data_dummy.dat','file');
        try;delete('mps_soft_data_dummy.dat');end
    end
    prior{ip}.MPS.soft_data_filename=['mps_soft_data_dummy.dat'];
end

%% RUN FORWARD
[m_propose{ip},prior{ip}.MPS]=mps_cpp(prior{ip}.ti,SIM,prior{ip}.MPS);

% In next iteration do NOT write the TI to disk!
prior{ip}.MPS.WriteTI=0;

%% Convert indexes to values
if isfield(prior{ip},'index_values');
    for i=1:length(prior{ip}.index_values)
        m_propose{ip}(find(m_propose{ip}==prior{ip}.index_values(i)))=prior{ip}.m_values(i);
    end
end


