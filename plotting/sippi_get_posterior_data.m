% sippi_get_posterior_data: load all data stored in mat-file
%
% Call: 
%
%  [data,prior,options,mcmc]=sippi_get_posterior_data(folder_name);
%  [data,prior,options,mcmc]=sippi_get_posterior_data(output_stuct);
%
function [data,prior,options,mcmc]=sippi_get_posterior_data(folder_name);

%options,prior,data,forward
%% LOAD THE CORRECT DATA
cwd=pwd;
if nargin==0
    % LOAD FROM MAT FILES
    [p,matfile]=fileparts(pwd);
    load(matfile);
elseif nargin==1;
    if isstruct(folder_name);
        folder_name=folder_name.txt;
    end
    fname=folder_name;
    cd(folder_name);
    load(folder_name);
end
if exist('C');
    mcmc=C{1}.mcmc;
end
cd(cwd);
    
