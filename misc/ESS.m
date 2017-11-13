% ESS: Effective Sample Size
%
% sample [nr,nm], nr:number of realizations, nm:number of model parameters
%
%  See also: multiESS
function [ess,nite_per_real]=ESS(sample,n_use,doPlot)


[nr,nm]=size(sample);

if nargin<2
    n_use = ceil(nr/4);
end
if n_use>nr
    n_use=nr-1;
    disp(sprintf('%s: n_use=%d',mfilename,n_use))
end
if nargin<3, doPlot=0;end


for im=1:nm;
    if nm>10
        progress_txt(im,nm)
    end
    s=sample(:,im);
    %ac2=acf(s,n_use);
    ac=autocorrelation(s);
    ac=ac(1:n_use); % do not use the whole time series
    
    nite_per_real(im) = (1+2*sum(ac)); 
    %sum(ac)
    %sum(ac2)
    ess(im)=nr/nite_per_real(im);    
    
    
    if doPlot==1
        plot(ac,'k-');hold on
        grid on
    end

        
end
if doPlot==1
    hold off
end