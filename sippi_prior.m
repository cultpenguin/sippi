
% sippi_prior A priori models for SIPPI 
%
% To generate a realization of the prior model defined by the prior structure use: 
%   [m_propose,prior]=sippi_prior(prior);
%
% To generate a realization of the prior model defined by the prior structure,
% in the vicinity of a current model (using sequential Gibbs sampling) use: 
%   [m_propose,prior]=sippi_prior(prior,m_current);
%
% The following types of a priori models can be used
%   % two point statistics bases
%   GAUSSIAN   [1D] : 1D generalized gaussian model
%   UNIFORM [1D-3D] : 1D-3D uncorrelated uniform distribution
%   CHOLESKY[1D-3D] : based on Cholesky decomposition
%   FFTMA   [1D-3D] : based on the FFT-MA method (Multivariate Gaussian) 
%   VISIM   [1D-3D] : based on Sequential Gaussian and Direct Sequential simulation
%   SISIM   [1D-3D] : based on Sequential indicator SIMULATION
%   % multiple point based statistics
%   SNESIM  [1D-3D] : based on a multiple point statistical model inferref from a training images. Relies in the SNESIM algorithm
%
%
%%%% SIMPLE EXAMPLE %%%
%
%% A simple 2D multivariate Gaissian based prior model based on the 
%% FFT-MA method, can be defined using 
%   im=1;
%   prior{im}.type='FFTMA';
%   prior{im}.name='A SIMPLE PRIOR';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%   prior{im}.m0=10;
%   prior{im}.Va='1 Sph(10)';
%   prior=sippi_prior_init(prior);
%% A realization from this prior model can be generated using
%   m=sippi_prior(prior);
%% This realization can now be plotted using
%   sippi_plot_prior(m,prior);
%% or
%   imagesc(prior{1}.x,prior{1}.y,m{1})
%
%%%% A PRIOR MODEL WITH SEVERAL 'TYPES OF A PRIORI MODEL'
%
%   im=1;
%   prior{im}.type='GAUSSIAN';
%   prior{im}.m0=100;
%   prior{im}.std=50;
%   prior{im}.norm=100;
%   im=im+1;
%   prior{im}.type='FFTMA';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%   prior{im}.m0=10;
%   prior{im}.Cm='1 Sph(10)';
%   im=im+1;
%   prior{im}.type='VISIM';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%   prior{im}.m0=10;
%   prior{im}.Cm='1 Sph(10)';
%   im=im+1;
%   prior{im}.type='SISIM';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%   prior{im}.m0=10;
%   prior{im}.Cm='1 Sph(10)';
%   im=im+1;
%   prior{im}.type='SNESIM';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%
%   sippi_plot_prior(prior);
%
%%% Sequential Gibbs sampling
%
%   All a priori model types can be perturbed, such that a new realization 
%   is generated in the vicinity of a current model. 
%   To do this Sequential Gibbs Sampling is used.
%   For more information, see <a href="matlab:web('http://dx.doi.org/10.1007/s10596-011-9271-1')">Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012. Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling. Computational Geosciences</a>.
%   The type of sequential Gibbs sampling can be controlled in the
%   'seq_gibbs' structures, e.g. prior{1}.seq_gibbs
%
%   im=1;
%   prior{im}.type='SNESIM';
%   prior{im}.x=[0:1:100];
%   prior{im}.y=[0:1:100];
%
%   [m,prior]=sippi_prior(prior);
%   prior{1}.seq_gibbs.step=1; % Large step--> independant realizations
%   prior{1}.seq_gibbs.step=.1; % Smaller step--> Dependant realizations
%   for i=1:30;
%      [m,prior]=sippi_prior(prior,m); % One iteration of Sequential Gibbs
%      sippi_plot_prior(prior,m);
%   end
%
% See also: sippi_prior_init, sippi_plot_prior, sippi_plot_prior_sample, sippi_prior_set_steplength.m
%
% TMH/2012
%

function [m_propose,prior]=sippi_prior(prior,m_current);



% Check for obsolete coding
for im=1:length(prior);
    if isfield(prior{im},'prior')
        if isfield(prior{im}.prior,'x');
            disp(sprintf('Using ''prior{im}.prior.x'' is now obsolete and should be ''prior{im}.x'' '))
            prior{im}.x=prior{im}.prior.x;
        end
        if isfield(prior{im}.prior,'y');
            disp(sprintf('Using ''prior{im}.prior.y'' is now obsolete and should be ''prior{im}.y'' '))
            prior{im}.x=prior{im}.prior.x;
        end
        if isfield(prior{im}.prior,'z');
            disp(sprintf('Using ''prior{im}.prior.z'' is now obsolete and should be ''prior{im}.z'' '))
            prior{im}.x=prior{im}.prior.x;
        end
        try
            disp(sprintf('Removing ''prior{im}.prior'''))
            prior{im}=rmfield(prior{im},'prior');
        end
    end
end
      

% Check for initialization
for im=1:length(prior);
    if ~isfield(prior{im},'init')
        prior=sippi_prior_init(prior);
    end
end


nm=length(prior);
if nargin>1
    m_propose=m_current;
end


%% SELECT WHICH MODEL PARAMETERS TO PERTURB
im_array=[];
for im=1:nm
    if isfield(prior{im},'perturb'),
        if prior{im}.perturb==1;
            im_array=[im_array im];
        end
    else
        im_array=[im_array im];
    end
    
end


if nargin==1;
    im_array=1:nm; % SAMPLE ALL MODEL PARAMETERS
end

if isempty(im_array)
    disp(sprintf('%s : no model parameters perturbed...',mfilename))
end

run_fftma=[];
% FIRST CHECK FOR ALL PRIOR TYPES, EXCEPT MASTER TYPES
for im=im_array;
    
    if isfield(prior{im},'Cm');
        if ~isfield(prior{im},'Va');
            prior{im}.Va=prior{im}.Cm;
        end
    end
      
    if  (strcmp(upper(prior{im}.type),'SISIM'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('sisim');
        end
        if isfield(prior{im},'marginal_prob');
            prior{im}.S.XML.parameters.Marginal_Probabilities.value=prior{im}.marginal_prob;
            % prior{im}.marginal_prob=ones(1,n_ind)/n_ind;
            n_ind=length(prior{im}.S.XML.parameters.Marginal_Probabilities.value);
            prior{im}.S.XML.parameters.Nb_Indicators.value=n_ind;
        end
        
        
        if isfield(prior{im},'Va');
            va_xml=sgems_variogram_xml(prior{im}.Va);
            prior{im}.S.XML.parameters.Variogram_Median_Ik=va_xml;
        end
        
        % NREAL = 1
        prior{im}.S.XML.parameters.Nb_Realizations.value=1;
        % SEED
        if isfield(prior{im},'seed');
            prior{im}.S.XML.parameters.Seed.value=prior{im}.seed;
        else
            prior{im}.S.XML.parameters.Seed.value=ceil(rand(1).*1e+6);
        end
        
        % REMOVE CONDITIONAL DATA.
        % FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
        if isfield(prior{im}.S,'f_obs')
            prior{im}.S=rmfield(prior{im}.S,'f_obs');
        end
        prior{im}.S.XML.parameters.Hard_Data_Grid.value='';
        prior{im}.S.XML.parameters.Hard_Data_Property.value='';
        
        %% SEQ GIBBS
        if nargin>1
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            [data,used]=sgems_set_resim_data(prior{im}.S,m_current{im},prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
            prior{im}.S=data; 
            prior{im}.seq_gibbs.used=used;
        end
        
        %% SISIM SIM
        prior{im}.S = sgems_grid(prior{im}.S);
        m_propose{im} = prior{im}.S.D';
        
    elseif  (strcmp(upper(prior{im}.type),'SNESIM')|isfield(prior{im},'S'))
        % SGEMS / SNESIM
        
        % REMOVE CONDITIONAL DATA.
        % FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
        if isfield(prior{im}.S,'f_obs')
            prior{im}.S=rmfield(prior{im}.S,'f_obs');
        end
        prior{im}.S.XML.parameters.Hard_Data.grid='';
        prior{im}.S.XML.parameters.Hard_Data.property='';
        
        %
        prior{im}.S.XML.parameters.Nb_Realizations.value=1;
        
        if isfield(prior{im},'seed');
            prior{im}.S.XML.parameters.Seed.value=prior{im}.seed;
        else
            prior{im}.S.XML.parameters.Seed.value=ceil(rand(1).*1e+6);
        end
        %disp(prior{im}.S.XML.parameters.Seed.value)
        % CHECK FOR SCALING AND ROTATION
          
        if isfield(prior{im},'rotation')
            prior{im}.S.XML.parameters.Use_Rotation.value=1;
            prior{im}.S.XML.parameters.Use_Global_Rotation.value=1;
            prior{im}.S.XML.parameters.Global_Angle.value=-1*prior{im}.rotation;
        end
        if isfield(prior{im},'scaling')
            if length(prior{im}.scaling)==1; aff=prior{im}.scaling.*[1 1 1]; end
            if length(prior{im}.scaling)==2; aff(3)=1; end
            prior{im}.S.XML.parameters.Use_Affinity.value=1;
            prior{im}.S.XML.parameters.Use_Global_Affinity.value=1;
            prior{im}.S.XML.parameters.Global_Affinity.value=aff;
        end
        
        if nargin>1
            % SEQUENTIAL GIBBS
            if isfield(prior{im},'index_values');
                m = zeros(size(m_current{im}))-1;
                for i=1:length(prior{im}.index_values)
                    try
                        m(find(m_current{im}==prior{im}.m_values(i)))=prior{im}.index_values(i);
                    end
                end
                m_current{im}=m;
            end
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            [data,used]=sgems_set_resim_data(prior{im}.S,m_current{im},prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
            prior{im}.S=data; 
            prior{im}.seq_gibbs.used=used;
        end
        
        prior{im}.S = sgems_grid(prior{im}.S);
        m_propose{im} = prior{im}.S.D';
        
        if isfield(prior{im},'index_values');
            for i=1:length(prior{im}.index_values)
                m_propose{im}(find(m_propose{im}==prior{im}.index_values(i)))=prior{im}.m_values(i);
            end
        end
        
        
    elseif (strcmp(upper(prior{im}.type),'VISIM')|isfield(prior{im},'V'))
        % VISIM PRIOR
        if isfield(prior{im},'Va');
            % update VISIM covariance settings
            if ~isstruct(prior{im}.Va);
                Va=deformat_variogram(prior{im}.Va);
            else
                Va=prior{im}.Va;
            end
            [prior{im}.V]=visim_set_variogram(prior{im}.V,prior{im}.Va);
            
        end
        
        if isfield(prior{im},'m0')
            prior{im}.V.gmean=prior{im}.m0; % set global mean
        end
        
        
        
        %% RANDOM SEED
        prior{im}.V.cond_sim=0;
        if isfield(prior{im},'seed');
            prior{im}.V.rseed=prior{im}.seed;
        else
            prior{im}.V.rseed=ceil(rand(1).*1e+6);
        end
        
        if isfield(prior{im},'d_target')
            f_cond=sprintf('d_target_%02d.eas',im);
            if ~exist(f_cond,'file');
                write_eas(f_cond,prior{im}.d_target);
            end
            prior{im}.V.refhist.fname=f_cond;
            prior{im}.V.ccdf=1;
            
            if isfield(prior{im},'min')
                prior{im}.V.tail.zmin=prior{im}.min;
            else
                prior{im}.V.tail.zmin=min(prior{im}.d_target);
            end
            if isfield(prior{im},'max')
                prior{im}.V.tail.zmin=prior{im}.max;
            else
                prior{im}.V.tail.zmax=max(prior{im}.d_target);
            end
                
            
        end
        
        %% SEQUENTIAL GIBBS
        if nargin>1
            % SEQUENTIAL GIBBS
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            %prior{im}.S=sgems_set_resim_data(prior{im}.S,m_current,prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
            [prior{im}.V, i_resim]=visim_set_resim_data(prior{im}.V,m_current{im},prior{im}.seq_gibbs.step,[],[],prior{im}.seq_gibbs.type);
            prior{im}.seq_gibbs.used=i_resim;
            if isempty(i_resim)
                prior{im}.V.cond_sim=0;
            else
                prior{im}.V.cond_sim=2;
            end
            
        end
        
        %% RUN VISIM
        prior{im}.V=visim(prior{im}.V);
        m_propose{im} = prior{im}.V.D';
        
        
    elseif (strcmp(lower(prior{im}.type),'gaussian'))
        %% 1D GENERALIZED GAUSSIAN
        
        if (isfield(prior{im},'d_target'))&(~isfield(prior{im},'o_nscore'))
            % UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
            d_min=min(prior{im}.d_target);
            d_max=max(prior{im}.d_target);
            [d_nscore,o_nscore]=nscore(prior{im}.d_target,1,1,d_min,d_max,0);
            prior{im}.o_nscore=o_nscore;
        end        
        
        if ~isfield(prior{im},'norm');
            prior{im}.norm=2;
        end
       
        if ~isfield(prior{im},'std');
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.std=(prior{im}.max-prior{im}.min)/2;
            else
                prior{im}.std=1;
            end
        end
        if ~isfield(prior{im},'m0');
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.m0=(prior{im}.max+prior{im}.min)/2;
            else
                prior{im}.m0=0;
            end
        end
        
        if prior{im}.norm==2
            if nargin>1
                if isfield(prior{im},'gauss_real')
                    gauss_real=prior{im}.gauss_real;
                else
                    gauss_real=(m_current{im}-prior{im}.m0)./prior{im}.std;
                end
                gauss_real_new=randn(1);
                step=pi/2*prior{im}.seq_gibbs.step;
                gauss_real=(cos(step)*gauss_real+sin(step)*gauss_real_new);
                
                %m_propose{im}=gauss_real.*prior{im}.std+prior{im}.m0;
            else
                gauss_real=randn(1);
            end
            prior{im}.gauss_real=gauss_real;
            if isfield(prior{im},'o_nscore')
                m_propose{im}=inscore(gauss_real,prior{im}.o_nscore);
            else
                m_propose{im}=gauss_real.*prior{im}.std+prior{im}.m0;
            end
            
        else
            if nargin>1
                if isfield(prior{im},'gauss_real');
                    gauss_real=prior{im}.gauss_real;
                else
                    gauss_real=randn(1);
                end
                gauss_real_new=randn(1);
                %t=sprintf('%5.3f %5.3f',gauss_real,gauss_real_new);
                step=pi/2*prior{im}.seq_gibbs.step;
                gauss_real=(cos(step)*gauss_real+sin(step)*gauss_real_new);
                %t=sprintf('%s %5.3f',t,gauss_real);
                %disp(t);
                
            else
                gauss_real=randn(1);
            end
            
            prior{im}.gauss_real=gauss_real;
            % transform gauss_real to generalized gauss
            g_cdf=normcdf(gauss_real,0,1);
            if ~isfield(prior{im},'ggauss_cdf');
                x=-5:.001:5;nx=length(x);
                pdf=generalized_gaussian(x,0,1,prior{im}.norm,0);
                prior{im}.ggauss_cdf.x=x;
                prior{im}.ggauss_cdf.cdf=cumsum(pdf)./sum(pdf);
            end
            
            % LOCATE VALUE IN GENERALIZED GAUSSIAN DIST.
            jj=find( prior{im}.ggauss_cdf.cdf>0 & prior{im}.ggauss_cdf.cdf<1);
            ggauss_real=interp1(prior{im}.ggauss_cdf.cdf(jj),prior{im}.ggauss_cdf.x(jj),g_cdf);
            if isnan(ggauss_real)
                ii=max(find(prior{im}.ggauss_cdf.cdf<g_cdf));
                if isempty(ii);ii=1;end
                ggauss_real=prior{im}.ggauss_cdf.x(ii);
            end
            
            
            m_propose{im}=prior{im}.m0+prior{im}.std*ggauss_real;
            
        end
       
        if ~isfield(prior{im},'round_ceil');
            prior{im}.round_ceil=0;
        end
        if prior{im}.round_ceil==1;
            m_propose{im}=ceil(m_propose{im});
        end
        
    elseif (strcmp(upper(prior{im}.type),'FFTMA'))        
        % THE FFTMA PRIOR IS HANDLED BELOW
        run_fftma=[run_fftma im];
    %% DEESSE (CLOSED SOURCE DIRECT SIMULATION
    elseif (strcmp(upper(prior{im}.type),'DS'))
        disp(sprintf('%s : ''%s'' type prior not yet implemented',mfilename,prior{im}.type));        
    else
        
        % check that a file exist that implements the prior type
        m_file=sprintf('sippi_prior_%s',lower(prior{im}.type));
        if exist(m_file,'file')
            p{1}=prior{im};
            if nargin==1
                [m_p,p]=feval(m_file,p);
            else
                m_c=m_current{im};
                [m_p,p]=feval(m_file,p,m_c);
            end
            m_propose{im}=m_p{1};
            prior{im}=p{1};
            
        else
            disp(sprintf('%s : ''%s'' type prior model not supported',mfilename,prior{im}.type));
        end
    end
        
end

%% CHECK IF WE NEED TO RUN FFTMA TYPE PRIOR BACUSE IT IS A MASTER 
run_fftma_as_master=[];
prior_master=[];

for im=im_array;
    
    if isfield(prior{im},'prior_master');
        try
        if (strcmp(upper(prior{prior{im}.prior_master}.type),'FFTMA'))
            run_fftma_as_master=[run_fftma_as_master prior{im}.prior_master];
        end
        prior_master=[prior_master prior{im}.prior_master];
        catch
            % perhaps prioir number prior{im}.prior_master does not exist
        end
    end
    prior_master=unique(prior_master);
    
end
run_fftma_as_master=unique(run_fftma_as_master);

%%
% WE NEED TO CHECK FOR FFTMA TYPE PRIOR SEPERATELY, AS IT CAN BE AFFECTED
% BY OTHER TYPES OF GAUSSIAN 1D TYPE PRIORS
im_fftma_array=unique([run_fftma_as_master,run_fftma]);
for im=im_fftma_array;
    if isempty(im); break;end % Needed as of Matlab R2013a
    if (strcmp(upper(prior{im}.type),'FFTMA'))
        %% FFTMA
        % update VA structure from range/ang/sill if possible
        if ~isstruct(prior{im}.Va); prior{im}.Va=deformat_variogram(prior{im}.Va);end
        
        % UPDATE COVARIANCE PARAMETERS IF THE HAVE BEEN DEFINED
        %[range,rot,sill,Va]=Va2RangeRot(prior{im}.Va);        
        %prior{im}.fftma_options.constant_C=1;
        
        % run thorugh the priors that have been updated
        im2_arr=im_array(1:(find(im_array==im)-1));
        im2_arr=1:length(prior);
        
        for j=im2_arr;
            
            update_master=0;
            % CHECK IF THIS PRIOR IS A MASTER
            if isfield(prior{j},'prior_master')
                if (prior{j}.prior_master==im);
                    % YES -> current im is master for current j'th prior
                    update_master=1;
                end
            end
            
            %if nargin==1, update_master=0;end
            if update_master==1
                
                % THIS IM IS A MASTER SO WE CAN UPDATE THE COVARIANCE MODEL
                % IF CHOSEN 
                if strcmp(prior{j}.name,'m0');
                    prior{im}.m0=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'sill');
                    prior{im}.Va.par1(1)=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_1');
                    range(1)=m_propose{j};
                    prior{im}.Va.par2(1)=range(1);
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_2');
                    range(2)=m_propose{j};
                    prior{im}.Va.par2(3)=range(2)/prior{im}.Va.par2(1);
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_3');
                    %range(3)=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'ang_1');
                    ang(1)=m_propose{j};
                    prior{im}.Va.par2(2)=ang(1);
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'nu');
                    nu=m_propose{j};
                    prior{im}.Va.nu=nu;
                    prior{im}.fftma_options.constant_C=0;
                end
            end
        end
        
        % set resim type and step length
        prior{im}.fftma_options.resim_type=prior{im}.seq_gibbs.type;
        prior{im}.fftma_options.lim=prior{im}.seq_gibbs.step;
        
        % IF CONSTANT_C=0, THEN REMOVE C and FFTC
        if prior{im}.fftma_options.constant_C==0;
            try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'fftC');end
            try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'C');end
        end
        
        %if prior{im}.perturb==0
        if isempty(run_fftma)|(run_fftma==0);
            % DO NOT PERTURB RANDOM NUMBERS UNLESS THE CURRENT FFT_MA TYPE
            % PRIOR IS ASKED TO BE PERTURBED
            % (ONLY COVARIANCE PROPERTIES ARE PERTURBED)
            prior{im}.fftma_options.lim=0;
        end
            
        
        %if (length(prior{im}.z)==1)&(length(prior{im}.y)==1)
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.Va,fftma_options);
        %elseif (length(prior{im}.z)==1)
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.Va,fftma_options);
        %else
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,fftma_options);
        %end
        [m_propose{im},z_rand,prior{im}.fftma_options]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,prior{im}.Va,prior{im}.fftma_options);
        prior{im}.fftma_options.z_rand=z_rand;
        
         if (isfield(prior{im},'d_target'))&(~isfield(prior{im},'o_nscore'))
              % UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
             d_min=min(prior{im}.d_target);
             d_max=max(prior{im}.d_target);
             [d_nscore,o_nscore]=nscore(prior{im}.d_target,1,1,d_min,d_max,0);
             prior{im}.o_nscore=o_nscore;
         end
        
        
        % PERFORM NORMAL SCORE OF NEEDED        
        if isfield(prior{im},'o_nscore');
            if ~isstruct(prior{im}.Va);
                prior{im}.Va=deformat_variogram(prior{im}.Va);
            end
            Va_par=prior{im}.Va;
            gvar=sum([Va_par.par1]);
            m_propose{im}=m_propose{im}./sqrt(gvar);
            m_propose{im}=inscore(m_propose{im},prior{im}.o_nscore);
        %else
        end
        
        % add mean model
        m_propose{im}=m_propose{im}+prior{im}.m0;
               
        prior{im}.m=m_propose{im};
    end
end


    
%% FIX EXTREME VALUES
for im=im_array;
    if isfield(prior{im},'min');
        try
        ii=find(m_propose{im}<prior{im}.min);
        catch
            keyboard
        end
        m_propose{im}(ii)=prior{im}.min;
        if nargin>2
            m_propose{im}(ii)=m_current{im}(ii);
        end
    end
    if isfield(prior{im},'max');
        ii=find(m_propose{im}>prior{im}.max);
        m_propose{im}(ii)=prior{im}.max;
        if nargin>2
            m_propose{im}(ii)=m_current{im}(ii);
        end
    end
end    

