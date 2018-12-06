% sippi_set_path Set paths for running sippi

function sippi_set_path();
%function [F,p]=sippi_set_path();
[p]=fileparts(which('sippi_set_path.m'));
if (isempty(p)|strcmp(p,'.'))
    p=pwd;
end


% adding toolboxes shipped with SIPPI
i=0;
i=i+1;F{i}=p;
i=i+1;F{i}=[p,filesep,'data',filesep,'crosshole'];
% i=i+1;F{i}=[p,filesep,'data',filesep,'ti'];
i=i+1;F{i}=[p,filesep,'plotting'];
i=i+1;F{i}=[p,filesep,'misc'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'fast_marching_kroon'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'fast_marching_kroon',filesep,'functions'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'fast_marching_kroon',filesep,'shortestpath'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'mpslib',filesep,'matlab'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'frequency_matching'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'fdem1d'];
i=i+1;F{i}=[p,filesep,'toolboxes',filesep,'fdem1d',filesep,'fdem1d'];

% ALSO ADD ALL DIRS IN 'toolboxes' folder not allready specified
toolboxes_dir=[p,filesep,'toolboxes'];
p_tb=dir(toolboxes_dir);
for i=1:length(p_tb)
    if (p_tb(i).isdir)&(~strcmp('..',p_tb(i).name))&(~strcmp('.',p_tb(i).name))
        dir_txt=[toolboxes_dir,filesep,p_tb(i).name];
        add_dir=1;
        for i_f=1:length(F);            
            if strcmp(F{i_f},dir_txt); add_dir=0; end
        end
        if add_dir==1           
            F{length(F)+1}=dir_txt;
            disp(sprintf('Adding folder %s',dir_txt))
        end
    end
end

% ACTUALLY ADD THE PATH
addpath(pwd);
for i=1:length(F);
	
    if exist(F{i},'dir')
        try
            addpath(F{i});
            stat='OK';
        catch
            stat='FAILED';
        end
        disp(sprintf('trying to add path %s [%s]',F{i},stat));
    else
        stat='NONEXISTENT';
        disp(sprintf('DIRECTORY DOES NOT EXIST : %s ',F{i}))
    end
    
end

try
    mgstat_set_path;
catch
    sippi_verbose(sprintf('%s : No path seems to be set of mGstat - SIPPI will probably fail',mfilename))
end


succ=savepath;
if succ==0
    sippi_verbose(sprintf('%s : saved path for later session',mfilename),0)
else
    sippi_verbose(sprintf('%s : COULD NOT SAVE PATH for later session',mfilename),0)
end

