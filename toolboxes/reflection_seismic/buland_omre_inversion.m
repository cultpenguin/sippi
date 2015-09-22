% buland_omre_inversion
%
% INPUT:
%    d: cell structure where each cell contains an NMO gather as a column, OR
%       on matric of size [3*ns,nx];
%
% 
% See also: buland_omre_setup
%
%
function  [M,V,vp_est,vs_est,rho_est,Cm_est]=buland_omre_inversion(d,log_vp0,log_vs0,log_rho0,angle,wavelet,Cd,Cm);

na=length(angle);

if iscell(d);
  nm_all=length(d{1});
  nx=length(d);
else
  [nm_all,nx]=size(d);
end
ns=nm_all/na;
disp(sprintf('%s: ns=%g',mfilename,ns))

% SETUP INVERSION MATRICES
[A,D,W]=buland_omre_setup(log_vp0,log_vs0,log_rho0,ns,angle,wavelet);
G=W*A*D;


m0=[ones(ns,1)*log_vp0;ones(ns,1)*log_vs0;ones(ns,1)*log_rho0];
  

% INVERT
for ix=1:nx
  
  progress_txt(ix,nx,'Buland and Omre linear NMO inversion')
  
  if iscell(d);
    d_obs=d{ix};
  else
    d_obs=d(:,ix);
  end
  
  % call solver each time
  %[m_est,Cm_est]=least_squares_inversion(G,Cm,Cd,m0,d_obs);
  
  % Much faster 
  if ix==1;
      S = Cd + G*Cm*G';
      K=(Cm*G')/S;
  end
  m_est  = m0 + K * (d_obs-G*m0);
  Cm_est = Cm - K * (G*Cm);
  
  M(:,ix)=m_est(:);
  V(:,ix)=diag(Cm_est);
 
end

vp_est=M(1:ns,:);
vs_est=M([1:ns]+ns,:);
rho_est=M([1:ns]+2*ns,:);


