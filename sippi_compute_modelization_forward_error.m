% sippi_compute_modelization_forward_error Computes an estimate of the modelization erro
%
% Computes and estimate of the Gaussian modelization error, N(dt,Ct)
% caused by the use of an imperfect forward kernel
%
% If called with only one output '[Ct]=sippi..]' then the Gaussian model is
% assumed by centered around 0, (dt{1}=0).
%
% Call
%   [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward_app,prior,data,N);
%
%
function [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward_app,prior,data,N);

if nargin<5,
    N=500;
end

N_app = length(forward_app);

try
    fname=sprintf('Ct_%s_N%d',forward_full.type,N);
catch
    fname=sprintf('Ct_N%d',N);
end
t_end_txt='';
t0=now;
%% Solve the forward problem N times for both choices of forward models
for i=1:N
    %progress_txt(i,N,'estimating forward response');
    progress_txt(i,N,sprintf('%s(%s)',mfilename,t_end_txt));
    if i==2,
        t0=now;
    end
    if i==1;
        % PRE ALLOCATE
        for j=1:length(data);
            try
                dd{j}=zeros(length(data{j}.i_use),N);
            catch
                dd{j}=zeros(length(data{j}.d_obs),N);
            end
        end    
    end
    
    
    m=sippi_prior(prior);
    
    [d_1,forward_full,prior,data]=sippi_forward(m,forward_full,prior,data);
    if N_app==1;
        [d_2{1},forward_app,prior,data]=sippi_forward(m,forward_app,prior,data);
    else
        for na=1:N_app;
            [d_2{na},forward_app{na},prior,data]=sippi_forward(m,forward_app{na},prior,data);
        end
    end
    
    
    
    for j=1:length(d_1);
        for na=1:N_app
            dd{j}(:,i,na)=d_1{j}-d_2{na}{j};
        end
    end
    
    if find(isnan(dd{j}(:,i,na))); keyboard;end
    
    if i>2
        t_end_txt = time_loop_end(t0,i-1,N-1);
    else
        t_end_txt = '';
    end
    
    
    txt_out=sprintf('%04d/%04d ESTIMATED TIME FOR COMPLETION : %s',i,N,t_end_txt);
    
    if (round(i/100)==(i/100))
        try
            if i==1;
                fid=fopen('ct_progress.txt','w');
            else
                fid=fopen('ct_progress.txt','a+');
            end
            fprintf(fid,'%s\n',txt_out);
            fclose(fid);
        end
        %    save(sprintf('%s_%d',fname,i));
    end
    
    
end

%% estimate modelization covariance models
for i=1:length(d_1);
    
    for na=1:N_app
        
        dt{na}{i}=nanmean(dd{i}(:,:,na)')';
        if nargout<2
            dt{na}{i}=dt{na}{i}.*0;
        end
        dd{i}(:,:,na)=dd{i}(:,:,na)-repmat(dt{na}{i},1,N);
        
        Ct{na}{i}=(1/N).*(dd{i}(:,:,na)*dd{i}(:,:,na)');
    end
end

%% RESHAPE OF N_app=1;
if N_app==1,
    Ct_org=Ct;
    dt_org=dt;
    clear Ct dt;
    for i=1:length(d_1);
        dt{i}=dt_org{1}{i};
        Ct{i}=Ct_org{1}{i};
    end
end


save(fname);






