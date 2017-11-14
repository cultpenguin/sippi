% ESS: Effective Sample Size
%
% call:
%   [ess,tau]=ESS(sample,n_use,doPlot,iLag)
% 
%   sample [nr,nm], nr:number of realizations, nm:number of model parameters
%
%  See also: multiESS
function [ess,tau]=ESS(sample,n_use,doPlot,iLag)


[nr,nm]=size(sample);

if nargin<2
    n_use = ceil(nr/4);
end
if isempty(n_use)
    n_use = ceil(nr/4);
end
if n_use>nr
    n_use=nr-1;
    disp(sprintf('%s: n_use=%d',mfilename,n_use))
end
if nargin<3, doPlot=0;end
if nargin<3, iLag=1;end


for im=1:nm;
    % if nm>10, progress_txt(im,nm), end
    s=sample(:,im);
    ac=autocorrelation(s);
    ac=ac(1:n_use); % do not use the whole time series
    
    tau(im) = (1+2*sum(ac)); 
    ess(im)=nr/tau(im);    
    
    if doPlot==1
        plot([1:length(ac)].*iLag,ac,'k-','LineWidth',.1);hold on
    end
        
end
if doPlot==1
    grid on
    xl=xlim;
    x_tau=linspace(xl(1),xl(2),51);
    h_tau=hist(iLag*tau,x_tau);
    h_tau=.2*h_tau/max(h_tau(:));
    bar(x_tau,h_tau,'r')
    %bar(x_tau,h_tau,'r-*')
    hold off
    xlabel('Lag')
    ylabel('Autocorrelation(lag)')
end