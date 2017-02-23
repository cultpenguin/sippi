% sippi_sequential_gibbs_resim: select model parameters for sequential
% gibbs resimulation
%
% Call
%
%  i_resim=sippi_sequential_gibbs_resim(prior,ip,type,step);
%
function i_resim=sippi_sequential_gibbs_resim(prior,ip);

if ~isfield(prior{ip},'init');
    prior=sippi_prior_init(prior);
end

if prior{ip}.seq_gibbs.type==1;
    % BOX TYPE
    
    
    i_resim=[];
    
    step = prior{ip}.seq_gibbs.step;
    
    if ~isfield(prior{ip}.seq_gibbs,'n_ite');
        prior{ip}.seq_gibbs.n_ite=1;
    end
    n_ite=prior{ip}.seq_gibbs.n_ite;
    
    if length(step)<2, step(2:3)=step(1);end
    if length(step)<3, step(3)=step(2);end
  
    ih=[];
    used=prior{ip}.xx.*0+1;
    for j=1:n_ite
        % 
        pos(1)=min(prior{ip}.x)+rand(1)*(max(prior{ip}.x)-min(prior{ip}.x));
        pos(2)=min(prior{ip}.y)+rand(1)*(max(prior{ip}.y)-min(prior{ip}.y));
        pos(3)=min(prior{ip}.z)+rand(1)*(max(prior{ip}.z)-min(prior{ip}.z));

        used(find( (abs(prior{ip}.xx-pos(1))<step(1)) & (abs(prior{ip}.yy-pos(2))<step(2)) ))=0;
    end
    i_resim=find(1-used);

    
    
    
    
    %sippi_verbose(sprintf('%s: BOX type resim not implemented',mfilename))
    
elseif prior{ip}.seq_gibbs.type==2; 
    % RANDOM MODEL PARAMETERS
    n_all=prod(prior{1}.dim);
    
    if prior{ip}.seq_gibbs.step<=1
        % use n_resim as a proportion of all random deviates
        n_resim=prior{ip}.seq_gibbs.step.*n_all;
    else
        n_resim=floor(prior{ip}.seq_gibbs.step); % such that (1.1)->(1.0)
    end
    n_resim=ceil(n_resim);
    n_resim = min([n_resim n_all]);
    i_resim=randomsample(n_all,n_resim);
    
    if isfield(prior{ip}.seq_gibbs,'template')
        % select model parameters in the vicininty of the 
        % allready chosen modelparameters using a template
        n_resim=length(i_resim);
        n_temp=size(prior{ip}.seq_gibbs.template,1);
        i_resim_t=zeros(1,n_temp*n_resim);
        
        k=0;
        usedim=[prior{ip}.dim(2) prior{ip}.dim(1) prior{ip}.dim(3)];
        for i=1:n_resim;
            if prior{ip}.ndim==1;
                sippi_verbose(sprintf('%s: teamplate note implemented for 1D prior',mfilename))
            else 
                ind=i_resim(i);
                %[ix,iy,iz]=ind2sub(prior{ip}.dim,ind);
                [ix,iy,iz]=ind2sub(usedim,ind);
                for j=1:n_temp
                    ix1=ix+prior{ip}.seq_gibbs.template(j,1);
                    iy1=iy+prior{ip}.seq_gibbs.template(j,2);
                    iz1=iz+prior{ip}.seq_gibbs.template(j,3);
                    try
                        ind_out = sub2ind(usedim,ix1,iy1,iz1);
                        k=k+1;
                        i_resim_t(k)=ind_out;
                    catch
                    end
                end
            end
        end
        
        % update i_resim
        i_resim=unique([i_resim i_resim_t(1:k)]);
        
    end
    
    
end



    








