% ESS: Effective Sample Size
%
% call:
%   [ess,tau]=ESS(sample,n_use,doPlot,iLag)
% 
%   sample [nr,nm], nr:number of realizations, nm:number of model parameters
%   n_use: The number of data point used to compute tau 
%          if n_use = 0, n_use os copmuted as 4 times the index of the
%          first negative value
%   doPlot [0/1]: Plot the autocorelation and Tau
%   iLag: The number of iteration between each realization, def=0;
%
%   if n_use=size(sample,1), then tau will tend be 0, and ess infinity
%
%  See also: multiESS
function [ess,tau,ess_all,tau_all]=ESS(sample,n_use,doPlot,iLag)


[nr,nm]=size(sample);

if nargin<2
    n_use = 0;
end
if n_use>nr
    n_use=nr-1;
    disp(sprintf('%s: n_use=%d',mfilename,n_use))
end
if nargin<3, doPlot=0;end
if nargin<4, iLag=1;end

%%
for im=1:nm;
    % if nm>10, progress_txt(im,nm), end
    s=sample(:,im);
    ac=autocorrelation(s);
    if n_use == 0;
        i_neg=find(ac<0);       
        if ~isempty(i_neg)
            n_use = i_neg(1)*4;
        else
            n_use = length(ac);
        end        
    end
    
    ac=ac(2:n_use); % do not use the whole time series, and do not use autocorrelation at lag=0;
    
    tau(im) = (1+2*sum(ac)); 
    ess(im)=nr/tau(im);    
    
    if doPlot==1
        plot([1:length(ac)].*iLag,ac,'k-','LineWidth',.1);hold on
    end
    
    if im==1;
        ac_all=ac;
    else
        ac_all = ac_all+ac;
    end
end
ac_all=ac_all/nm;
tau_all = (1+2*sum(ac_all));
ess_all=nr/tau_all;
 

plot([1:length(ac)].*iLag,ac_all,'g-','LineWidth',2);hold on
if doPlot==1
    grid on
    xl=xlim;
    x_tau=linspace(xl(1),xl(2),51);
    h_tau=hist(iLag*tau,x_tau);
    h_tau=.2*h_tau/max(h_tau(:));
    bar(x_tau,h_tau,'r')
    %bar(x_tau,h_tau,'r-*')
    plot([1,1].*tau_all*iLag,ylim,'g:')
    hold off
    xlabel('Lag')
    ylabel('Autocorrelation(lag)')
end


tau=tau*iLag;
tau_all=tau_all*iLag;

