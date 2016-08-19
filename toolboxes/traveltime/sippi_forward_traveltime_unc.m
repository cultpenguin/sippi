% sippi_forward_traveltime_unc: as sippi_forward_traveltime while updateing uncertainty in data
%
% Performs exactly as '' expect that the uncorrelated uncertainty 
% on data is allowed to changed
%
% To set the noise accroding to 1D prior distritbution, define a prior
% structrue with name 'd_std'.
% Then the standard devaition of the uncertainty will be see using this
% prior. 
%     
%   prior{1}.name=d_std';
%
% If no such prior is set, 'sippi_forward_traveltime' will be run.
%
% Call :
%   [d,forward,prior,data]=sippi_forward_traveltime_unc(m,forward,prior,data)
%
% See also: sippi_forward_traveltime
function [d,forward,prior,data]=sippi_forward_traveltime_unc(m,forward,prior,data,id_array,im)
if nargin<6,
    try
        im=length(prior);
    catch
        im=1;
    end
end
if nargin<5,
    try
        id_array=1:length(data);
    catch
        id_array=1;
    end
end
if nargin<2, forward.null=[]; end
if nargin<4, data{1}.null='';end

if nargin==2;
    [d,forward]=sippi_forward_traveltime(m,forward);
elseif nargin==3
    [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior);
elseif nargin==4
    [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data);
end

if nargin>=4
    for im=1:length(prior);
        if strcmp(lower(prior{im}.name),'d_std');
            d_std=m{im};
            sippi_verbose(sprintf('%s: updating noise on data to d_std=%6.2g',mfilename,d_std),2);
            for id=1:length(data)
                data{id}.d_std=d_std;
                if isfield(data{id},'Cd');
                    data{id}=rmfield(data{id},'Cd');
                end
                data{id}.d_std=d_std;
                
                data{id}.recomputeCD=1;
                data{id}.full_likelihood=1;
            end
        end
    end
end    
    
