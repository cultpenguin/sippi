function [data,prior,options,mcmc]=sippi_get_posterior_data(options);
%options,prior,data,forward
%% LOAD THE CORRECT DATA
cwd=pwd;
if nargin==0
    % LOAD FROM MAT FILES
    [p,matfile]=fileparts(pwd);
    load(matfile);
elseif nargin==1;
    if isstruct(options),
    else
        fname=options;
        cd(fname);
        load(fname);
    end
end
if exist('C');
    mcmc=C{1}.mcmc;
end
cd(cwd);
    
