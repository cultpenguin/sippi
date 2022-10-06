% sippi_compute_modelization_forward_error Computes an estimate of the modelization erro
%
% Computes and estimate of the Gaussian modelization error, N(dt,Ct)
% caused by the use of an imperfect forward kernel
%
% If called with only one output '[Ct]=sippi..]' then the Gaussian model is
% assumed by centered around 0, (dt{1}=0).
%
% Call
%   [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward_app,prior,N,d);
%
%
% For details see:
%  Hansen, T.M., Cordua, K. S., Jacobsen, B. J., and Mosegaard, K. (2014)
%  Accounting for imperfect forward modeling in geophysical inverse problems - exemplified for cross hole tomography.
%  Geophsyics, 79(3) H1-H21, 2014. doi:10.1190/geo2013-0215.1
%
function [Ct,dt,dd,d_full,d_app,o_nscore,dd_org]=sippi_compute_modelization_forward_error(forward_full,forward_app,prior,N,d,data,useNormalScore);


N_app = length(forward_app);

if N_app==1;
    if isstruct(forward_app);
        f_temp=forward_app;
        clear forward_app;
        forward_app{1}=f_temp;
    end
end
if nargin<4,N=500;end
if nargin<5
    % precalcualte once (to initialize forward structure) 
    m=sippi_prior(prior);
    %if N_app==1;
    %    [d,forward_app]=sippi_forward(m,forward_app,prior);
    %else
    for i=1:N_app
        [d{i},forward_app{i}]=sippi_forward(m,forward_app{i},prior);
    end
    %end
else
    if iscell(d);
        d_temp=d;
        clear d;
        d{1}=d_temp;
    end
end
if nargin<6
    data{1}.d_obs=[]
    data{1}.d_std=zeros(1,length(d{1}));
end
if nargin<7
    useNormalScore = 0;
    o_nscore{1}=[];;
end

try
    fname=sprintf('Ct_%s_N%d',forward_full.type,N);
catch
    fname=sprintf('Ct_N%d',N);
end
t_end_txt='';
time_per_step=1e-19;
t0=now;
%% Solve the forward problem N times for both choices of forward models
% PRE ALLOCATE

ND=length(d{1}{1});
if nargout>3; d_full{1}=zeros(ND,N); end
if nargout>4;
    for j=1:N_app;
        d_app{j}=zeros(ND,N);
    end
end
for j=1:N_app;
    dd{j}=zeros(ND,N);
end

n_retry=0;
i=0;
while i<N
    i=i+1;

    try
        
        % display progress
        txt=sprintf('%s(end->%s)',mfilename,t_end_txt);
        if (time_per_step>.5);
            disp(sprintf('%03d/%03d: %s',i,N,txt));
        else
            if (i/25)==round(i/25);
                progress_txt(i,N,txt);
            end
        end
        
        % generate realization from prior
        m=sippi_prior(prior);
        
        % solve FULL forward problem
        [d_1,forward_full,prior]=sippi_forward(m,forward_full,prior);
        
        % SIMULATE AN ERROR / ONLY FOR TESTING
        %if rand(1)>.8;
        %    d_1{1}=d_1{1}(1:10);
        %end
        
        
        % solve APPROXIMATE forward problem(s)
        for na=1:N_app;
            [d_2{na},forward_app{na},prior]=sippi_forward(m,forward_app{na},prior);
        end
        
        % compute and store realization of modeling error
        if nargout>3
            d_full{1}(:,i)=d_1{1};
        end
        for na=1:N_app
            dd{na}(:,i)=d_1{1}-d_2{na}{1};
            if nargout>4
                d_app{na}(:,i)=d_2{na}{1};
            end
        end
        
        % compute time to end of loop
        if i>1
            t_end_txt = time_loop_end(t0,i-1,N-1);
        end
        
        % write estimated time for end of simulat
        txt_out=sprintf('%04d/%04d ESTIMATED TIME FOR COMPLETION : %s',i,N,t_end_txt);
        if (round(i/100)==(i/100))
            % solve FULL forward problem
            try
                if (i/100)==1;
                    fid=fopen('ct_progress.txt','w');
                else
                    fid=fopen('ct_progress.txt','a+');
                end
                fprintf(fid,'%s\n',txt_out);
                fclose(fid);
            end
        end
        
        n_retry=0;
    
    catch
        
        n_retry_max=10;
        n_retry=n_retry+1;
        
        if n_retry>=n_retry_max;
            disp(sprintf('%s: Something went wrong - %d TIMES - STOPPING ',mfilename,n_retry));
            return
        end
        
        disp(sprintf('%s: Something went wrong - trying again (%02d/%d) ',mfilename,n_retry,n_retry_max));
        %i_old=i;
        i=max([1 i-1]);
        %disp(sprintf('%d -> %d',i_old,i))
    end
    
    time_per_step=((now-t0)*3600*24)/i;
    
end

%% estimate modelization covariance models

for na=1:N_app
    
    if useNormalScore==1;
        
        % simulate measurement noise
        dd_meas{na} = 0.*dd{na};
        if isfield(data{1},'Cd')
            dd_meas{na} = gaussian_simulation_cholesky(data{1}.d0, data{1}.Cd,N)
        elseif isfield(data{1},'d_var')
            for i=1:N
                dd_meas{na}(:,i) = randn(size(data{1}.d_var)).*sqrt(data{1}.d_var);
            end
        elseif isfield(data{1},'d_std')
            for i=1:N
                dd_meas{na}(:,i) = randn(size(data{1}.d_std)).*data{1}.d_std;
            end
        end

        nd=size(dd{na},1);
        % Convert sample of modeling PLUS measurement error into NS space
        for id=1:nd
            %[d_nscore{na}(id,:),o_nscore{na}{id}]=nscore(dd{na}(id,:));
            dd_org{na}(id,:)=dd{na}(id,:) + dd_meas{na}(id,:);
            [d_nscore{na}(id,:),o_nscore{na}{id}]=nscore(dd_org{na}(id,:));
        end
        % test forward
        [d_nscore_forward{na}]=nscore_mul(dd_org{na},o_nscore{na});
        % d_nscore_forward{na}  == d_nscore{na}
        dd{na}=d_nscore{na};
        
    end

    dt{na}=nanmean(dd{na}')';
    % optionally force dt=0;
    if nargout<2
        dt{na}=dt{na}.*0;
    end
    dd_ct{na}=dd{na}-repmat(dt{na},1,N);
    Ct{na}=(1/N).*(dd_ct{na}*dd_ct{na}');
    
end


save(fname);






