% sippi_get_resim_data: Get conditional data for resimulation
%
% d_cond=sippi_get_resim_data(m_current,prior,ip);
% 
% c_cond [n_cond,4]: col1: x, col2: y, col4: z, col4: d
%
% See also sippi_prior, sippi_sequential_gibbs_resim
%
function d_cond=sippi_get_resim_data(m_current,prior,ip);
if nargin<3
    ip=1;
end

% get index of data to resimulate
i_resim=sippi_sequential_gibbs_resim(prior,ip);
%disp(sprintf('N_Resim=%d/%d\n',length(i_resim),prod(size(m_current{ip}))));
% get index of hard conditional data
ih=setxor(1:(prod(prior{ip}.dim)),i_resim);

d_cond=[prior{ip}.xx(ih(:)) prior{ip}.yy(ih(:)) prior{ip}.zz(ih(:)) m_current{ip}(ih(:))];
