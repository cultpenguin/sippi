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
% See also: pg_plot
%
function [pg_d,pg_limits,pg_arr]=pg_transform(m,pg_map,pg_limits);

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
    n1=pg_size(1);
    w1=diff(pg_limits)/n1;
    pg_arr{1}=pg_limits(1)+w1/2+w1.*(0:1:[n1-1])  ;
    
    di_1=interp1(pg_arr{1},1:length(pg_arr{1}),m{1},'nearest','extrap');
    
    % use index to get value for pg_map. must be a faster way
    pg_d=di_1.*NaN;
    
    [ny,nx,nz]=size(di_1);
    for iz=1:nz
        for ix=1:nx
            for iy=1:ny
                pg_d(iy,ix,iz)=pg_map(di_1(iy,ix,iz));
            end
        end
    end
    
elseif pg_dim==2
    %% PLURIGAUSSIAN USING TWO GAUSSIANS
    % set the arrays
    n1=pg_size(2);
    n2=pg_size(1);
    
    w1=diff(pg_limits)/n1;
    pg_arr{1}=pg_limits(1)+w1/2+w1.*(0:1:[n1-1])  ;
    w2=diff(pg_limits)/n2;
    pg_arr{2}=pg_limits(1)+w2/2+w2.*(0:1:[n2-1])  ;
    
    % get the index
    di_1=interp1(pg_arr{1},1:length(pg_arr{1}),m{1},'nearest','extrap');
    di_2=interp1(pg_arr{2},1:length(pg_arr{2}),m{2},'nearest','extrap');    
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
    n1=pg_size(2);
    n2=pg_size(1);
    n3=pg_size(3)
    
    
    w1=diff(pg_limits)/n1;
    pg_arr{1}=pg_limits(1)+w1/2+w1.*(0:1:[n1-1])  ;
    w2=diff(pg_limits)/n2;
    pg_arr{2}=pg_limits(1)+w2/2+w2.*(0:1:[n2-1])  ;
    w3=diff(pg_limits)/n3;
    pg_arr{3}=pg_limits(1)+w3/2+w3.*(0:1:[n3-1])  ;
    
    di_1=interp1(pg_arr{1},1:length(pg_arr{1}),m{1},'nearest','extrap');
    di_2=interp1(pg_arr{2},1:length(pg_arr{2}),m{2},'nearest','extrap');
    di_3=interp1(pg_arr{3},1:length(pg_arr{3}),m{3},'nearest','extrap');
    
    [ny,nx,nz]=size(di_1);
    for iz=1:nz
        for ix=1:nx
            for iy=1:ny                
                pg_d(iy,ix,iz)=pg_map(di_2(iy,ix,iz),di_1(iy,ix,iz),di_3(iy,ix,iz));
            end
        end
    end
    
    
end
    