% sippi_least_squares Least squares type inversion for SIPPI
%
% Call :
%    [options,data,prior,forward,m_reals,m_est,Cm_est]=sippi_least_squares(data,prior,forward,options);
%
%   options.lsq.type    : LSQ type to use ('lsq' (classical linear leqast squares) is the default)
%   options.lsq.n_reals : Number of realizations to generate
%   options.lsq.plot : [0/1] show figures or not def->0.
%   options.lsq.save_data : [0/1] save realizations to output folder. def->1.
%
%
% TMH/01/2017
%
% See also sippi_rejection, sippi_metropolis
%



%              'error_sim', simulation through error simulation
%              'visim', simulation through SGSIM of DSSIM
%
function [m_est_out,options,data,prior,forward]=sippi_tikhonov(data,prior,forward,options);

id=1;
im=1;

m_reals=[];
m_est=[];
Cm_est=[];
options.lsq.null=''; % make sutre options.lsq exists


% save data?
if ~isfield(options.lsq,'save_data')
    options.lsq.save_data=1;
end

% interactive?
if ~isfield(options.lsq,'interactive')
    options.lsq.interactive=0;
end

% plot?
if ~isfield(options.lsq,'plot')
    options.lsq.plot=0;
end


%% CHOOSE NAME
if ~isfield(options,'txt')
    options.txt=mfilename;
end

try
    options.txt=sprintf('%s_%s_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt,options.lsq.type);
catch
    options.txt=sprintf('%s_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);
end
sippi_verbose(sprintf('%s: output folder: %s ',mfilename,options.txt),1)


%% MODEL COVARINCE
if ~isfield(options.lsq,'Cm');
    prior=sippi_prior_init(prior);
    if isfield(prior{im},'Cmat');
        options.lsq.Cm=prior{im}.Cmat;
    else
        prior=sippi_prior_init(prior);
        options.lsq.Cm=precal_cov([prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],[prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],prior{im}.Va);
    end
end
if ~isfield(options.lsq,'Cm');
    sippi_verbose(sprintf('%s: Could not model covariance Cm. Please use a Gaussian prior or set prior{%d}.Cmat',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: Model covariance, Cm, set in options.lsq.Cm',mfilename),1)
end

%% DATA COVARINCE
if isfield(data{id},'CD');
    options.lsq.Cd=data{id}.CD;
else
    % solve the forwrad problem and make use of data{1}.CD if it exists
    m=sippi_prior(prior);
    [d,forward,prior,data]=sippi_forward(m,forward,prior,data);
    [logL,L,data]=sippi_likelihood(d,data,id);
    try
        options.lsq.Cd=data{id}.CD;
    end
end

if ~isfield(options.lsq,'Cd');
% No correlated noise is set...
if isfield(data{id},'d_std');
    if length(data{id}.d_std)==1;
        options.lsq.Cd=eye(length(data{id}.d_obs)).*data{id}.d_std.^2;
    else
        options.lsq.Cd=diag(data{1}.d_std.^2);
    end
end
if isfield(data{id},'d_var');
    if length(data{id}.d_var)==1;
        options.lsq.Cd=eye(length(data{id}.d_obs)).*data{id}.d_var;
    else
        options.lsq.Cd=diag(data{1}.d_var);
    end
end
end


if ~isfield(options.lsq,'Cd');
    sippi_verbose(sprintf('%s: Could not data covariance Cd. Please use a Gaussian noise model in data{%d}',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: Data covariance, Cd, set in options.lsq.Cd',mfilename),1)
end


%% CHECK FOR FORWARD OPERATOR
if ~isfield(forward,'G');
    % assume the forward operator is output in forward.G if sippi_forward
    % is run
    try
        m=sippi_prior(prior);
        if isfield(forward,'forward_function');
            [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
        else
            [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im);
        end
    end
    if ~isfield(forward,'G');
        sippi_verbose(sprintf('%s : No forward operator G found in forward',mfilename),0)
    end
end
try
    options.lsq.G=forward.G;
end
if ~isfield(options.lsq,'G');
    sippi_verbose(sprintf('%s: linear forward operator G is not set. Please set in in forward.G',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: linear forward operator G set in options.lsq.G',mfilename),1)
end


%% M
if ~isfield(prior{im},'m0');
    prior{im}.m0=0;
end
if length(prior{im}.m0)==1;
    nm=size(options.lsq.Cm,1);
    options.lsq.m0=ones(nm,1).*prior{im}.m0;
else
    options.lsq.m0=prior{im}.m0(:);
end
sippi_verbose(sprintf('%s: setting options.lsq.m0=prior{%d}.m0',mfilename,im),1)

%% D
if ~isfield(data{id},'d0');
    options.lsq.d0=data{id}.d_obs.*0;
    sippi_verbose(sprintf('%s: setting options.lsq.d0=0;',mfilename),1)
end
options.lsq.d_obs=data{id}.d_obs;
sippi_verbose(sprintf('%s: setting options.lsq.d_obs=data{%d}.d_obs',mfilename,id),1)


%% THIKONOV
if ~isfield(options.lsq,'tikhonov');options.lsq.tikhonov=4;end
options.txt=[options.txt,'_',num2str(options.lsq.tikhonov)];
    

x = prior{1}.x;nx=length(x);
y = prior{1}.y;ny=length(y);

I = eye(size(options.lsq.Cm,1));

e_N = 21;
e_min = -3;
e_max = 3;
e_array = logspace(e_min,e_max, e_N);
for i=1:length(e_array);
    progress_txt(i,e_N,mfilename)
    e = e_array(i);
    
    if options.lsq.tikhonov==1
        % classical
        options.lsq.m0=0;
        m_est{i} = inv(options.lsq.G'*options.lsq.G + e.*I)*options.lsq.G'*(options.lsq.d_obs(:));
    elseif options.lsq.tikhonov==2
        % general, m0
        m_est{i} = options.lsq.m0 + inv(options.lsq.G'*options.lsq.G + e^2.*I)*options.lsq.G'*(options.lsq.d_obs(:)-options.lsq.G*options.lsq.m0);
    elseif options.lsq.tikhonov==3
        % general, m0, Cd, weigthed
        m_est{i} = options.lsq.m0 + inv(options.lsq.G'*options.lsq.Cd*options.lsq.G + e^2.*I)*options.lsq.G'*options.lsq.Cd*(options.lsq.d_obs(:)-options.lsq.G*options.lsq.m0);
    elseif options.lsq.tikhonov==4
        % general, m0, Cd, Cm
        m_est{i} = options.lsq.m0 + inv(options.lsq.G'*options.lsq.Cd*options.lsq.G + e^2.*inv(options.lsq.Cm))*options.lsq.G'*options.lsq.Cd*(options.lsq.d_obs(:)-options.lsq.G*options.lsq.m0);   
    end
    
    m_norm(i) = norm(m_est{i} - options.lsq.m0);
    d_norm(i) = norm(options.lsq.d_obs - options.lsq.G*m_est{i});
    
end

%% PLOT
figure(3);
plot(m_norm,d_norm,'k-*')
xlabel('||m-m0||')
ylabel('||d-Gm||')

% FIND OPTIMAL E
if options.lsq.interactive == 1;
    figure(3);
    pos=ginput(1);
    dd=abs(d_norm-pos(2));
    ic=find(dd==min(dd));
    e_opt = e_array(ic)
    
else
    % find max curvuature
    c=abs(curvuature_2d(m_norm,d_norm));
    ic=find(c==max(c));ic=ic(1);
    e_opt = e_array(ic)
    %figure(1);clf
    %semilogx(e_array,c);
end
% find best e using noise level
norm_noise = sqrt(length(data{1}.d_obs))*data{1}.d_std(1);
dd_norm = abs(norm_noise-d_norm);
ic_noise=find(dd_norm==min(dd_norm));ic_noise=ic_noise(1);
e_opt_noise =  e_array(ic_noise)


figure(3);
hold on
plot(m_norm(ic),d_norm(ic),'k.','MarkerSize',42)
text(m_norm(ic)*1.15,d_norm(ic)*1.15,sprintf('\\alpha = %3.1e',e_array(ic)))
hold off
print_mul(sprintf('%s_Lcurve',options.txt))

hold on
norm_noise = sqrt(length(data{1}.d_obs))*data{1}.d_std(1);
xl=xlim;
plot(xl,[1 1].*norm_noise,'b--');
hold off
print_mul(sprintf('%s_Lcurve_N',options.txt))



%%

m_est_out{1} = m_est{ic};

%%
figure(2);subfigure(1,1,1);clf;
N_plot=5;
i_plot=ceil(linspace(1,e_N,N_plot));
for j=1:N_plot;
    i=i_plot(j);
    subplot(1,N_plot+2,j)
    imagesc(x,y,reshape(m_est{i},ny,nx));
    axis image
    try;caxis(prior{1}.cax);end
    title(sprintf('\\alpha=%3.1e',e_array(i)))
end
figure(2);
subplot(1,N_plot+2,N_plot+1)
imagesc(x,y,reshape(m_est{ic},ny,nx));
axis image
try;caxis(prior{1}.cax);end
title(sprintf('Optimal L \n e=%3.1e',e_opt))

subplot(1,N_plot+2,N_plot+2)
imagesc(x,y,reshape(m_est{ic_noise},ny,nx));
axis image
try;caxis(prior{1}.cax);end
title(sprintf('Optimal (noise)\n e=%3.1e',e_opt_noise))


print_mul(sprintf('%s_mest',options.txt))


% ALSO PLOT THE MODEL THAT MATCH THE NOISE������

%%



%% LOG LOG
%% find max curvuature
c2=curvuature_2d(log10(m_norm),log10(d_norm));
ic2=find(c2==max(c2));ic2=ic2(1);
e_opt2 = e_array(ic2);
figure(4);
semilogx(m_norm,d_norm,'-*')
hold on
plot(m_norm(ic2),d_norm(ic2),'k.','MarkerSize',42)
text(m_norm(ic2)*1.05,d_norm(ic2)*1.01,sprintf('\\alpha = %3.1e',e_array(ic2)))
hold off
xlabel('||m-m0||')
ylabel('||d-Gm||')
print_mul(sprintf('%s_logLcurve',options.txt))    

save T
