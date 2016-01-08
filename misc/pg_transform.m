% pg_transform: plurigaussian transformation
%
% Call: 
%    [pg_d]=pg_transform(m,pg_map,pg_limits);
%      m: realizations of gaussian distribtion(s)
%      pg_map: Map defining pluri-Gaussian transformation 
%      pg_limits: Limit values for pg_map (def:[-3 3])
%
%      pd_d: Pluri-Gaussian transformed data
% 
% See also: sippi_prior_plurigaussian, prior_reals_plurigaussian
%
% Example
%
function [pg_d]=pg_transform(m,pg_map,pg_limits);

%% check input data
if nargin<2
    pg_map=[0 1];
end
if nargin<3
    pg_limits=[-3 3];
end

% m must be a cell
if ~iscell(m);
    d_tmp=m;
    clear m
    m{1}=d_tmp;
end
if ~iscell(m)
    disp(sprintf('%s: first input must be a cell when pg_dim>1',mfilename));
end


% find the 'dimension' of the transformation
pg_dim=ndims(pg_map);
pg_size=size(pg_map);

% check/fix for 1D
if find(pg_size==1)
    pg_map=pg_map(:);
    pg_size=max(pg_size);
    pg_dim=1;
end
    


%% 

if pg_dim==1;
    %% TRUNCATED GAUSSIAN
    pg_arr_1=linspace(pg_limits(1),pg_limits(2),pg_size(1));    
    di_1=interp1(pg_arr_1,1:length(pg_arr_1),m{1},'nearest','extrap');
    
    % use index to get value for pg_map. must be a faster way
    pg_d=di_1.*NaN;
    
    [ny,nx,nz]=size(di_1);
    for iz=1:nz
        for ix=1:nx
            for iy=1:ny
                pg_d(iy,ix,iz)=pg_map(di_1(iy,ix));
            end
        end
    end
    
elseif pg_dim==2
    %% PLURIGAUSSIAN USING TWO GAUSSIANS
    % set the arrays
    pg_arr_1=linspace(pg_limits(1),pg_limits(2),pg_size(2));
    pg_arr_2=linspace(pg_limits(1),pg_limits(2),pg_size(1));

    % get the index
    di_1=interp1(pg_arr_1,1:length(pg_arr_1),m{1},'nearest','extrap');
    di_2=interp1(pg_arr_2,1:length(pg_arr_2),m{2},'nearest','extrap');
    
    [ny,nx,nz]=size(di_1);
    for iz=1:nz
        for ix=1:nx
            for iy=1:ny                
                pg_d(iy,ix,iz)=pg_map(di_2(iy,ix,iz),di_1(iy,ix,iz));
            end
        end
    end
    
    
elseif pg_dim==3
    %% PLURIGAUSSIAN USING THREE GAUSSIANS
    pg_arr_1=linspace(pg_limits(1),pg_limits(2),pg_size(2));
    pg_arr_2=linspace(pg_limits(1),pg_limits(2),pg_size(1));
    pg_arr_3=linspace(pg_limits(1),pg_limits(2),pg_size(3));
    
    di_1=interp1(pg_arr_1,1:length(pg_arr_1),m{1},'nearest','extrap');
    di_2=interp1(pg_arr_2,1:length(pg_arr_2),m{2},'nearest','extrap');
    di_3=interp1(pg_arr_3,1:length(pg_arr_3),m{3},'nearest','extrap');
    
    [ny,nx,nz]=size(di_1);
    for iz=1:nz
        for ix=1:nx
            for iy=1:ny                
                pg_d(iy,ix,iz)=pg_map(di_2(iy,ix,iz),di_1(iy,ix,iz),di_2(iy,ix,iz));
            end
        end
    end
    
    
end
    