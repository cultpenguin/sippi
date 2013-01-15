% sippi_set_path Set paths for running sippi

function sippi_set_path();
[p]=fileparts(which('sippi_set_path.m'));
if isempty(p)
    p=pwd;
end

i=0;
i=i+1;F{i}='';
i=i+1;F{i}=['data',filesep,'crosshole',filesep];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon',filesep];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon',filesep,'functions',filesep];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon',filesep,'shortestpath',filesep];
i=i+1;F{i}=['toolboxes',filesep,'traveltime',filesep];
i=i+1;F{i}=['toolboxes',filesep,'lomgrav',filesep];
%i=i+1;F{i}=['toolboxes',filesep,'reflection',filesep];
%i=i+1;F{i}=['toolboxes',filesep,'gpr_full_waveform',filesep];

%if exist([p,filesep,F{i}],'dir')
%    % add path to full waveform GPR forward
%    i=i+1;F{i}=['toolboxes',filesep,'gpr_full_waveform',filesep,'forward_simulation',filesep];
%    for  j=1:18;
%        i=i+1;F{i}=['toolboxes',filesep,'gpr_full_waveform',filesep,'forward_simulation',filesep,sprintf('Core%d',j),filesep];
%    end
%end
i=i+1;F{i}=['toolboxes',filesep,'mGstat',filesep];



addpath(pwd);
for i=1:length(F);
    path_str=sprintf('%s%s%s',p,filesep,F{i});
    
    if exist(path_str,'dir')
        try
            addpath(path_str);
            stat='OK';
        catch
            stat='FAILED';
        end
        disp(sprintf('trying to add path %s [%s]',path_str,stat));
    else
        stat='NONEXISTENT';
        disp(sprintf('DIRECTORY DOES NOT EXIST : %s ',path_str))
    end
    
end

try
    mgstat_set_path;
catch
    disp(sprintf('%s : No path seems to be set of mGstat - SIPPI will probably fail',mfilename))
end