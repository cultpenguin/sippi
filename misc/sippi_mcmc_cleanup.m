% sippi_mcmc_cleanup: remove output dir of sippi_metropolis
%
% sippi_mcmc_cleanup(options.out);
% sippi_mcmc_cleanup(options.out);
%
%
function sippi_mcmc_cleanup(fname)

if isstruct(fname)
    fname=fname.txt;      
end


if ~exist(fname,'dir')
    sippi_verbose(sprintf('%s: folder ''%s'' does not exist',mfilename,fname),0);
    return
end
sippi_verbose(sprintf('%s: trying to delete folder ''%s''',mfilename,fname),0);

f=dir(fname);
for i_file=1:length(f);
    if f(i_file).isdir==0
        fdel = fullfile(f(i_file).folder,f(i_file).name);
        sippi_verbose(sprintf('%s: trying to delete file ''%s''',mfilename,fdel),0);
        delete(fdel)
    end
end

[status, message, messageid]  = rmdir(fname,'s');





