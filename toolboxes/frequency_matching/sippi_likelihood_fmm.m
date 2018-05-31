% sippi_likelihood_multinomial: Compute likelihood using the multinomial function
%
% Call:
%   [logL,L,data]=sippi_likelihood_fmm(d,data);
%
% Input parameter:
%    d{id}; Frequency distribution to be evaluated. 
%    data{id}.d_obs; Observed frequency distribution.
%    data{id}.nprior; Prior frequency distribution.
%
% This likelihood uses the multinomial distribution to calculate the
% probability of some (e.g. pattern) frequency distribution given some
% observed frequency distribution and a prior frequency distribution.
%
% See also sippi_likelihood
%
% For more details on the theory see:
% Cordua et al., 2015. Improving the Pattern Reproducibility of Multiple-Point-Based Prior
% Models Using Frequency Matching. Mathematical Geosciences. April 2015,
% Volume 47, Issue 3, pp 317?343. DOI: 10.1007/s11004-014-9531-4.1111


function [logL,logL_all,data]=sippi_likelihood_fmm(d,data,id_array)
if  nargin<3
    LL=length(d);
    id_array=1:LL;
else
    LL=length(length(id_array));
end



logL_all=zeros(1,LL);
for id=id_array
    
    [~,X]=size(d{id});
    if X>1
        disp('d{id} has to be a column vector')
        return
    end
    [~,X]=size(data{id}.d_obs);
    if X>1
        disp('data{id}.d_obs has to be a column vector')
        return
    end
    %logL_all(id)=0;

    if ~isfield(data{id},'nprior');
        data{id}.nprior=0;
    end
   
    logL_all(id)=multinomial(d{id},data{id}.d_obs,data{id}.nprior);
            
end
logL=sum(logL_all);

