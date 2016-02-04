% setup_Cm_corr: setup correlated covariance model structure according to 
%                Buland and Omre (2003)
%
%
% Example:
%
%  mu_omre = log([3000 2250 3000]);
%  var_omre = [0.0074 0.00240 .0074 ];
%  cc=[1 0.8 0.7;0.8 1 0.5;0.7 0.5 1];
%  Va=sprintf('0.001 Nug(0) + 0.999 Gau(20)');
%  pos=[1:1:80]';
%  Cm0=precal_cov(pos,pos,Va);
%  [m0,Cmat]=setup_Cm_corr(Cm0,mu_omre,var_omre,cc);
%
%
% See also; sippi_prior_cholesky
%
function [m,Cmat]=setup_Cm_corr(Cm0,m0,var0,cc);
  N=length(var0);
          
  nCm=size(Cm0,1);
  m=zeros(N*nCm,1);

  Cmat=[];
  for i=1:N;
      m([1:nCm]+(i-1)*nCm)=m0(i);
      % loop over row
      Cmat_row=[];
      for j=1:N;
          fak_i = cc(i,j)*sqrt(var0(i))*sqrt(var0(j));
          if j==1;
              Cmat_row=Cm0*fak_i;
          else
              Cmat_row=[Cmat_row Cm0*fak_i];
          end
      end
      if i==1;
          Cmat=Cmat_row;
      else
          Cmat=[Cmat;Cmat_row];
      end
  end
        
