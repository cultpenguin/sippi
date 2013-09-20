% sippi_set_path Set paths for running sippi

function sippi_set_path();
[p]=fileparts(which('sippi_set_path.m'));
if isempty(p)
    p=pwd;
end


i=0;
i=i+1;F{i}='';
i=i+1;F{i}=['data',filesep,'crosshole'];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon'];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon',filesep,'functions'];
i=i+1;F{i}=['toolboxes',filesep,'fast_marching_kroon',filesep,'shortestpath'];
%i=i+1;F{i}=['toolboxes',filesep,'traveltime'];
%i=i+1;F{i}=['toolboxes',filesep,'lomgrav'];
%i=i+1;F{i}=['toolboxes',filesep,'mGstat'];
%if exist([p,filesep,F{i}],'dir')
%    % add path to full waveform GPR forward
%    i=i+1;F{i}=['toolboxes',filesep,'gpr_full_waveform',filesep,'forward_simulation'];
%    for  j=1:18;
%        i=i+1;F{i}=['toolboxes',filesep,'gpr_full_waveform',filesep,'forward_simulation',filesep,sprintf('Core%d',j)];
%    end
%end


% ALSO ADD ALL DIRS IN 'toolboxes' folder not allready specified
toolboxes_dir=[p,filesep,'toolboxes'];
p_tb=dir(toolboxes_dir)
for i=1:length(p_tb)
    if (p_tb(i).isdir)&(~strcmp('..',p_tb(i).name))&(~strcmp('.',p_tb(i).name))
        dir_txt=['toolboxes',filesep,p_tb(i).name];
        add_dir=1;
        for i_f=1:length(F);            
            if strcmp(F{i_f},dir_txt); add_dir=0; end
        end
        if add_dir==1           
            F{length(F)+1}=dir_txt;
            disp(sprintf('AAADDING %s',dir_txt))
        end
    end
end

% ACTUALLY ADD THE PATH
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