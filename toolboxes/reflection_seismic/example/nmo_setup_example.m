%% nmo_setup_prior;
rng('default');
rng(1);

% PRIOR TYPE
if ~exist('ptype','var');  ptype=[2]; end
% SIGNAL TO NOISE
if ~exist('SN','var');  SN=1; end

% Set nmo gather center locations
dx=10;
if ~exist('nx','var');nx=11;end
x=[0:1:(nx-1)].*dx;

% Set time
dt=0.001;
t=2.1:dt:2.3;
nt=length(t);

%% SETUP PRIOR STRUCTURE

% Covariance model properties
h_x=80;
ang=90;
h_t=0.01;
h_aniso=h_t/80;

% mean and standard devaition of elastic parameters
vp_m0=3000;
vp_std=260;

vs_m0=3000;
vs_std=190;

rho_m0=2250;
rho_std=75;

% three types of prior model
% ptype=1; ->> FFTMA vp,vs,rho
% ptype=2; ->> FFTMA log(vp),log(vs), and log(rho)
% ptype=3; ->> CHOLESKY [log(vp);log(vs);log(rho)]


if ptype==1
  % VP VS RHO
  ip=1;
  prior{ip}.name='vp';
  prior{ip}.type='fftma';
  prior{ip}.x=x;
  prior{ip}.y=t;
  prior{ip}.m0=vp_m0;
  prior{ip}.Cm=sprintf('%g Exp(%f,90,%3.3g)',vp_std^2,h_x,h_t/h_x);
  
  ip=ip+1;
  prior{ip}.name='vs';
  prior{ip}.type='fftma';
  prior{ip}.x=x;
  prior{ip}.y=t;
  prior{ip}.m0=vs_m0;
  prior{ip}.Cm=sprintf('%g Exp(%f,90,%3.3g)',vs_std^2,h_x,h_t/h_x);
    
  ip=ip+1;
  prior{ip}.name='rho';
  prior{ip}.type='fftma';
  prior{ip}.x=x;
  prior{ip}.y=t;
  prior{ip}.m0=rho_m0;
  prior{ip}.Cm=sprintf('%g Exp(%f,90,%3.3g)',rho_std^2,h_x,h_t/h_x);
  
  forward.log=0;
  
elseif ptype==2
  
  
  % BULAND and OMRE type log parameterization
  mu_omre = log([vp_m0 vs_m0 rho_m0 ]);
  var_omre = [0.0074 0.00240 .0074 ];
  r_omre=20; % ms
  
  ip=0;
  ip=ip+1;
  prior{ip}.name='vp';
  prior{ip}.type='fftma';
  prior{ip}.x=x;
  prior{ip}.y=t;
  %prior{ip}.y=2100:1:2105;ny=length(prior{ip}.y);
  prior{ip}.m0=(mu_omre(ip));
  prior{ip}.Cm=sprintf('%g Gau(%f,90,%3.3g)',var_omre(ip),h_x,h_t/h_x);
  
  ip=ip+1;
  prior{ip}=prior{ip-1};
  prior{ip}.name='vs';
  prior{ip}.m0=(mu_omre(ip));
  prior{ip}.Cm=sprintf('%g Gau(%f,90,%3.3g)',var_omre(ip),h_x,h_t/h_x);
 
  ip=ip+1;
  prior{ip}=prior{ip-1};
  prior{ip}.name='rho';
  prior{ip}.m0=(mu_omre(ip));
  prior{ip}.Cm=sprintf('%g Gau(%f,90,%3.3g)',var_omre(ip),h_x,h_t/h_x);
   
  forward.log=1;
  
end

%% SETUP WAVELETS
forward.angle=[0:5:45]; % angle for each trace in NMO gather
forward.freq=linspace(60,30,length(forward.angle)); % frequency of each gather

%% SETUP FORWARD MODEL
forward.forward_function='sippi_forward_reflection_nmo';
forward.t=t;
forward.type='zoeppritz';

%% REFERENCE MODEL
[m_ref,prior]=sippi_prior(prior);
figure(1);
for i=1:3;
  subplot(1,3,i);
  imagesc(m_ref{i});title(prior{i}.name);colorbar
end

%% REFERENCE DATA
forward.type='zoeppritz';
[d_ref,forward]=sippi_forward(m_ref,forward,prior);

use_corr_noise=1;
S=std(d_ref{1}(:));
  
if use_corr_noise==1;
  % Correlated data noise
  h_a=5;
  h_t=0.01;
  d_var=(S/SN).^2;
  Va=sprintf('.00001 Nug(0) + %f Gau(%g,90,%f)',d_var,h_a,h_t/h_a);
  da=1;
  na=length(forward.angle);
  Cd=precal_cov_2d(na,nt,da,dt,Va);
end

for id=1:length(x);
  data{id}.d_ref=d_ref{id};
  
  if use_corr_noise==1;
    % corr noise
    data{id}.Cd=Cd;
    data{id}.d_noise=gaussian_simulation_cholesky(0,Cd,1);
  else   
    % uncorr noise
    data{id}.d_std=S/SN;
    data{id}.d_noise=randn(size(d_ref{id})).*data{id}.d_std;
  end
  data{id}.d_obs=data{id}.d_ref+data{id}.d_noise;
end

sippi_plot_data_reflection_nmo(d_ref,data,1,prior);

save(sprintf('nmo_reference_data_type%d_SN%g_nx%02d',ptype,SN,nx));
